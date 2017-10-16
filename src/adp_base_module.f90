! ------------------------------------------------------------------------
! Copyright (C) 2017 Nobuya Sato, Hiori Kino, and Takashi Miyake
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
! ------------------------------------------------------------------------

! Angular-dependent potential (ADP)
! Y. Mishin, M. J. Mehl, and D. A. Papaconstantopoulos, Acta Mater. 53, 4029 (2005).
module adp_base_module

  use, non_intrinsic :: periodic_lattice_module, only: &
      & cell_list_type, &
      & periodic_lattice_type

  implicit none

  private

  integer, parameter, public :: kr = kind(0d0)

  type, abstract, public :: adp_base_type

    private

    type(periodic_lattice_type) :: lattice
    integer, allocatable :: z(:)
    real(kr), allocatable :: r_frac(:, :)
    type(cell_list_type), allocatable :: cell_of_replica(:, :)
    real(kr), allocatable :: rhobar(:)
    real(kr), allocatable :: mu(:, :)
    real(kr), allocatable :: lambda(:, :, :)

  contains

    private

    procedure, non_overridable, public :: set_system
    procedure, non_overridable, public :: energy
    procedure, non_overridable, public :: force
    procedure, non_overridable, public :: stress

    procedure(cutoff_interface), deferred :: cutoff ! cutoff of interaction
    procedure(phi_interface), deferred :: phi
    procedure(phi_interface), deferred :: phi_derivative
    procedure(f_interface), deferred :: f
    procedure(f_interface), deferred :: f_derivative
    procedure(rho_interface), deferred :: rho
    procedure(rho_interface), deferred :: rho_derivative
    procedure(u_interface), deferred :: u
    procedure(u_interface), deferred :: u_derivative
    procedure(w_interface), deferred :: w
    procedure(w_interface), deferred :: w_derivative

    procedure, non_overridable :: phi_force
    procedure, non_overridable :: psi_force

    procedure, nopass, non_overridable :: nu

  end type adp_base_type

  abstract interface

    pure real(kr) function cutoff_interface(this)
      import :: adp_base_type
      import :: kr
      class(adp_base_type), intent(in) :: this
    end function cutoff_interface

    pure real(kr) function phi_interface(this, z1, z2, r)
      import :: adp_base_type
      import :: kr
      class(adp_base_type), intent(in) :: this
      integer, intent(in) :: z1
      integer, intent(in) :: z2
      real(kr), intent(in) :: r
    end function phi_interface

    pure real(kr) function f_interface(this, z, rhobar)
      import :: adp_base_type
      import :: kr
      class(adp_base_type), intent(in) :: this
      integer, intent(in) :: z
      real(kr), intent(in) :: rhobar
    end function f_interface

    pure real(kr) function rho_interface(this, z, r)
      import :: adp_base_type
      import :: kr
      class(adp_base_type), intent(in) :: this
      integer, intent(in) :: z
      real(kr), intent(in) :: r
    end function rho_interface

    pure real(kr) function u_interface(this, z1, z2, r)
      import :: adp_base_type
      import :: kr
      class(adp_base_type), intent(in) :: this
      integer, intent(in) :: z1
      integer, intent(in) :: z2
      real(kr), intent(in) :: r
    end function u_interface

    pure real(kr) function w_interface(this, z1, z2, r)
      import :: adp_base_type
      import :: kr
      class(adp_base_type), intent(in) :: this
      integer, intent(in) :: z1
      integer, intent(in) :: z2
      real(kr), intent(in) :: r
    end function w_interface

  end interface

contains

  subroutine set_system(this, a, z, r_frac)

    use, intrinsic :: iso_fortran_env, only: error_unit
    use, non_intrinsic :: periodic_lattice_module, only: krpl => kr

    implicit none

    class(adp_base_type), intent(inout) :: this
    real(kr), intent(in) :: a(3, 3) ! direct lattice vectors as column vectors
    integer, intent(in) :: z(:) ! atomic numbers
    real(kr), intent(in) :: r_frac(:, :) ! fractional atomic coords.

    real(kr), allocatable :: rhobar(:)
    real(kr), allocatable :: mu(:, :)
    real(kr), allocatable :: lambda(:, :, :)
    integer :: num_atom
    real(kr) :: rj_frac_reduced(3)
    real(kr) :: rj_frac(3)
    real(kr) :: rij(3)
    real(kr) :: dij
    integer :: i
    integer :: j
    integer :: jj

    if (size(r_frac, dim=1) /= 3 .or. size(r_frac, dim=2) /= size(z)) then
      write(error_unit, "(a, ': ', a, ': ', a)") "adp_base_type", "set_system", "shape of 'r_frac' is not 3 x 'size(z)'"
      stop 1
    end if

    this%lattice%direct_lattice = real(a, krpl)
    this%z = z
    this%r_frac = r_frac

    num_atom = size(z)

    if (allocated(this%cell_of_replica)) then
      deallocate(this%cell_of_replica)
    end if
    allocate(this%cell_of_replica(num_atom, num_atom))

    !$omp parallel do private(j)
    do i = 1, num_atom
      do j = 1, num_atom
        this%cell_of_replica(j, i) = this%lattice%get_cell_of_replica( &
            & real(r_frac(:, j), krpl), &
            & real(r_frac(:, i), krpl), &
            & real(this%cutoff(), krpl))
      end do
    end do
    !$omp end parallel do

    allocate(rhobar(num_atom))
    allocate(mu(3, num_atom))
    allocate(lambda(3, 3, num_atom))

    rhobar = 0._kr
    mu = 0._kr
    lambda = 0._kr
    !$omp parallel do reduction(+: rhobar, mu, lambda), private(rj_frac_reduced, rj_frac, rij, dij, j, jj)
    do i = 1, num_atom
      do j = 1, num_atom

        rj_frac_reduced = modulo(r_frac(:, j), 1._kr)

        do jj = 1, size(this%cell_of_replica(j, i)%list, dim=2)

          if (j == i .and. all(this%cell_of_replica(j, i)%list(:, jj) == floor(r_frac(:, i)))) then
            cycle
          end if

          rj_frac = this%cell_of_replica(j, i)%list(:, jj) + rj_frac_reduced
          rij = matmul(a, rj_frac - r_frac(:, i))
          dij = norm(rij)

          ! Eq. (2)
          rhobar(i) = rhobar(i) + this%rho(z(j), dij)

          ! Eq. (3)
          mu(:, i) = mu(:, i) + this%u(z(i), z(j), dij) * rij

          ! Eq. (4)
          lambda(:, :, i) = lambda(:, :, i) + this%w(z(i), z(j), dij) * direct_product(rij, rij)

        end do

      end do
    end do
    !$omp end parallel do

    this%rhobar = rhobar
    this%mu = mu
    this%lambda = lambda

    deallocate(lambda)
    deallocate(mu)
    deallocate(rhobar)

  end subroutine set_system

  ! Eq. (1)
  real(kr) function energy(this)

    use, intrinsic :: iso_fortran_env, only: error_unit

    implicit none

    class(adp_base_type), intent(in) :: this

    real(kr) :: e
    integer :: num_atom
    real(kr) :: rj_frac_reduced(3)
    real(kr) :: rj_frac(3)
    real(kr) :: rij(3)
    real(kr) :: dij
    integer :: i
    integer :: j
    integer :: jj

    if (.not.allocated(this%z)) then
      write(error_unit, "(a, ': ', a, ': ', a)") "adp_base_type", "energy", "system is not set"
      stop 1
    end if

    num_atom = size(this%z)

    e = 0._kr

    ! 1st term
    !$omp parallel do reduction(+: e), private(rj_frac_reduced, rj_frac, rij, dij, j, jj)
    do i = 1, num_atom
      do j = 1, num_atom

        rj_frac_reduced = modulo(this%r_frac(:, j), 1._kr)

        ! sum over periodic replicas of the j-th atom
        do jj = 1, size(this%cell_of_replica(j, i)%list, dim=2)

          if (j == i .and. all(this%cell_of_replica(j, i)%list(:, jj) == floor(this%r_frac(:, i)))) then
            cycle
          end if

          rj_frac = this%cell_of_replica(j, i)%list(:, jj) + rj_frac_reduced
          rij = matmul(this%lattice%direct_lattice, rj_frac - this%r_frac(:, i))
          dij = norm(rij)

          e = e + 0.5_kr * this%phi(this%z(i), this%z(j), dij)

        end do

      end do
    end do
    !$omp end parallel do

    ! 2nd term
    !$omp parallel do reduction(+: e)
    do i = 1, num_atom
      e = e + this%f(this%z(i), this%rhobar(i))
    end do
    !$omp end parallel do

    ! 3rd term
    e = e + 0.5_kr * sum(this%mu ** 2)

    ! 4th term
    e = e + 0.5_kr * sum(this%lambda ** 2)

    ! 5th term
    !$omp parallel do reduction(+: e)
    do i = 1, num_atom
      e = e - this%nu(this%lambda(:, :, i)) ** 2 / 6._kr
    end do
    !$omp end parallel do

    energy = e

  end function energy

  ! Appendix
  function force(this)

    use, intrinsic :: iso_fortran_env, only: error_unit

    implicit none

    class(adp_base_type), intent(in) :: this
    real(kr), allocatable :: force(:, :)

    real(kr), allocatable :: f(:, :)
    integer :: num_atom
    real(kr) :: rj_frac_reduced(3)
    real(kr) :: rj_frac(3)
    real(kr) :: rij(3)
    integer :: i
    integer :: j
    integer :: jj

    if (.not.allocated(this%z)) then
      write(error_unit, "(a, ': ', a, ': ', a)") "adp_base_type", "force", "system is not set"
      stop 1
    end if

    num_atom = size(this%z)

    allocate(f(3, num_atom))

    f = 0._kr
    !$omp parallel do reduction(+: f), private(rj_frac_reduced, rj_frac, rij, j, jj)
    do i = 1, num_atom
      do j = 1, num_atom

        rj_frac_reduced = modulo(this%r_frac(:, j), 1._kr)

        ! sum over periodic replicas of the j-th atom
        do jj = 1, size(this%cell_of_replica(j, i)%list, dim=2)

          if (j == i .and. all(this%cell_of_replica(j, i)%list(:, jj) == floor(this%r_frac(:, i)))) then
            cycle
          end if

          rj_frac = this%cell_of_replica(j, i)%list(:, jj) + rj_frac_reduced
          rij = matmul(this%lattice%direct_lattice, rj_frac - this%r_frac(:, i))

          ! 1st term of the 1st eq. in Appendix
          f(:, i) = f(:, i) + this%phi_force(this%z(i), this%z(j), this%rhobar(i), this%rhobar(j), rij)

          ! 2nd term of the 1st eq. in Appendix
          f(:, i) = f(:, i) + this%psi_force( &
              & this%z(i), this%z(j), &
              & this%mu(:, i), this%mu(:, j), &
              & this%lambda(:, :, i), this%lambda(:, :, j), &
              & rij)

        end do

      end do
    end do
    !$omp end parallel do

    force = f

    deallocate(f)

  end function force

  function stress(this)

    use, intrinsic :: iso_fortran_env, only: error_unit

    implicit none

    class(adp_base_type), intent(in) :: this
    real(kr) :: stress(3, 3)

    real(kr) :: s(3, 3)
    real(kr) :: f(3)
    integer :: num_atom
    real(kr) :: rj_frac_reduced(3)
    real(kr) :: rj_frac(3)
    real(kr) :: rij(3)
    integer :: i
    integer :: j
    integer :: jj

    if (.not.allocated(this%z)) then
      write(error_unit, "(a, ': ', a, ': ', a)") "adp_base_type", "stress", "system is not set"
      stop 1
    end if

    num_atom = size(this%z)

    s = 0._kr
    !$omp parallel do reduction(+: s), private(rj_frac_reduced, rj_frac, rij, f, j, jj)
    do i = 1, num_atom
      do j = 1, num_atom

        rj_frac_reduced = modulo(this%r_frac(:, j), 1._kr)

        ! sum over periodic replicas of the j-th atom
        do jj = 1, size(this%cell_of_replica(j, i)%list, dim=2)

          if (j == i .and. all(this%cell_of_replica(j, i)%list(:, jj) == floor(this%r_frac(:, i)))) then
            cycle
          end if

          rj_frac = this%cell_of_replica(j, i)%list(:, jj) + rj_frac_reduced
          rij = matmul(this%lattice%direct_lattice, rj_frac - this%r_frac(:, i))

          ! 1st term of the 1st eq. in Appendix
          f = this%phi_force(this%z(i), this%z(j), this%rhobar(i), this%rhobar(j), rij)

          ! 2nd term of the 1st eq. in Appendix
          f = f + this%psi_force( &
              & this%z(i), this%z(j), &
              & this%mu(:, i), this%mu(:, j), &
              & this%lambda(:, :, i), this%lambda(:, :, j), &
              & rij)

          s = s + 0.5_kr * direct_product(f, -rij) / this%lattice%volume()

        end do

      end do
    end do
    !$omp end parallel do

    stress = s

  end function stress

  ! 2nd eq. in Appendix
  pure function phi_force(this, zi, zj, rhobar_i, rhobar_j, rij)

    implicit none

    class(adp_base_type), intent(in) :: this
    integer, intent(in) :: zi
    integer, intent(in) :: zj
    real(kr), intent(in) :: rhobar_i
    real(kr), intent(in) :: rhobar_j
    real(kr), intent(in) :: rij(3)
    real(kr) :: phi_force(3)

    real(kr) :: dij

    dij = norm(rij)

    ! 1st term originating from the 1st term of Eq. (1)
    phi_force = this%phi_derivative(zi, zj, dij) * rij / dij

    ! 2nd term originating from the 2nd term of Eq. (1)
    phi_force = phi_force &
        & + (this%f_derivative(zi, rhobar_i) * this%rho_derivative(zj, dij) &
        &    + this%f_derivative(zj, rhobar_j) * this%rho_derivative(zi, dij)) * rij / dij

  end function phi_force

  ! 3rd eq. in Appendix
  pure function psi_force(this, zi, zj, mu_i, mu_j, lambda_i, lambda_j, rij)

    implicit none

    class(adp_base_type), intent(in) :: this
    integer, intent(in) :: zi
    integer, intent(in) :: zj
    real(kr), intent(in) :: mu_i(3)
    real(kr), intent(in) :: mu_j(3)
    real(kr), intent(in) :: lambda_i(3, 3)
    real(kr), intent(in) :: lambda_j(3, 3)
    real(kr), intent(in) :: rij(3)
    real(kr) :: psi_force(3)

    real(kr) :: dij
    real(kr) :: w
    real(kr) :: w_derivative

    dij = norm(rij)
    w = this%w(zi, zj, dij)
    w_derivative = this%w_derivative(zi, zj, dij)

    ! 1st and 2nd terms originating from the 3rd term of Eq. (1)
    psi_force = (mu_i - mu_j) * this%u(zi, zj, dij)
    psi_force = psi_force + dot_product((mu_i - mu_j), rij) * this%u_derivative(zi, zj, dij) * rij / dij

    ! 3rd and 4th terms originating from the 4th term of Eq. (1)
    psi_force = psi_force + 2._kr * matmul(rij, lambda_i + lambda_j) * w
    psi_force = psi_force + dot_product(rij, matmul(lambda_i + lambda_j, rij)) * w_derivative * rij / dij

    ! 5th term originating from the 5th term of Eq. (1)
    ! The correct sign is negative. See the footnote 1 in Y. Mishin and A. Lozovoi, Acta Mater. 54, 5013 (2006).
    psi_force = psi_force - (this%nu(lambda_i) + this%nu(lambda_j)) * (w_derivative * dij + 2._kr * w) * rij / 3._kr

  end function psi_force

  ! Eq. (5)
  pure real(kr) function nu(lambda)

    implicit none

    real(kr), intent(in) :: lambda(3, 3)

    nu = lambda(1, 1) + lambda(2, 2) + lambda(3, 3)

  end function nu

  pure real(kr) function norm(v)

    implicit none

    real(kr), intent(in) :: v(:)

    norm = sqrt(sum(v ** 2))

  end function norm

  pure function direct_product(v1, v2)

    implicit none

    real(kr), intent(in) :: v1(:)
    real(kr), intent(in) :: v2(:)
    real(kr), allocatable :: direct_product(:, :)

    integer :: j

    allocate(direct_product(size(v1), size(v2)))

    do j = 1, size(v2)
      direct_product(:, j) = v1 * v2(j)
    end do

  end function direct_product

end module adp_base_module
