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

! Lennard-Jones potential
! N. Karasawa and W. A. Goddard III, J. Phys. Chem. 93, 7320 (1989).
module lj_module

  use, non_intrinsic :: periodic_lattice_module, only: &
      & cell_list_type, &
      & periodic_lattice_type

  implicit none

  private

  integer, parameter, public :: kr = kind(0d0)

  real(kr), parameter :: pi = 4._kr * atan(1._kr)
  real(kr), parameter :: sqrt_pi = sqrt(pi)

  type, public :: lj_parameter_type
    private
    integer, allocatable, public :: z(:)
    real(kr), allocatable, public :: epsilon(:, :)
    real(kr), allocatable, public :: sigma(:, :)
  end type lj_parameter_type

  type, public :: lj_type
    private
    type(lj_parameter_type), public :: param
    type(periodic_lattice_type) :: lattice
    integer, allocatable :: z(:)
    real(kr), allocatable :: r_frac(:, :)
    real(kr), allocatable :: a(:, :)
    real(kr), allocatable :: b(:, :)
    real(kr) :: energy_error_tolerance
    real(kr) :: eta
    type(cell_list_type), allocatable :: cell_real(:, :)
    type(cell_list_type) :: cell_reciprocal
  contains
    private
    procedure, public :: set_system
    procedure, public :: energy
    procedure, public :: force
    procedure, public :: stress
    procedure, public :: energy_error_real
    procedure, public :: energy_error_reciprocal
  end type lj_type
  interface lj_type
    module procedure construct_lj_type
  end interface lj_type

  public :: make_parameter_lorentz_berthelot
  public :: make_parameter_kong
  public :: make_parameter_waldman_hagler

contains

  function construct_lj_type(param, tol) result(lj)

    use, intrinsic :: iso_fortran_env, only: error_unit

    implicit none

    type(lj_parameter_type), intent(in) :: param
    real(kr), intent(in), optional :: tol
    type(lj_type) :: lj

    if (.not. allocated(param%z)) then
      write(error_unit, "(a, ': ', a)") "lj_type", "'param%z' is not allocated"
      stop 1
    end if
    if (.not. allocated(param%epsilon)) then
      write(error_unit, "(a, ': ', a)") "lj_type", "'param%epsilon' is not allocated"
      stop 1
    end if
    if (.not. allocated(param%sigma)) then
      write(error_unit, "(a, ': ', a)") "lj_type", "'param%sigma' is not allocated"
      stop 1
    end if

    if (any(shape(param%epsilon) /= size(param%z))) then
      write(error_unit, "(a, ': ', a)") "lj_type", "sizes of 'param%epsilon' does not match the size of 'param%z'"
      stop 1
    end if
    if (any(shape(param%sigma) /= size(param%z))) then
      write(error_unit, "(a, ': ', a)") "lj_type", "sizes of 'param%sigma' does not match the size of 'param%z'"
      stop 1
    end if

    if (any(param%epsilon /= transpose(param%epsilon))) then
      write(error_unit, "(a, ': ', a)") "lj_type", "'param%epsilon' is not symmetric"
      stop 1
    end if
    if (any(param%sigma /= transpose(param%sigma))) then
      write(error_unit, "(a, ': ', a)") "lj_type", "'param%sigma' is not symmetric"
      stop 1
    end if

    if (any(param%epsilon < 0._kr)) then
      write(error_unit, "(a, ': ', a)") "lj_type", "'param%epsilon' is not positive"
      stop 1
    end if
    if (any(param%sigma < 0._kr)) then
      write(error_unit, "(a, ': ', a)") "lj_type", "'param%sigma' is not positive"
      stop 1
    end if

    lj%param = param

    if (present(tol)) then
      if (tol < 0._kr) then
        write(error_unit, "(a, ': ', a)") "lj_type", "'tol' is not positive"
        stop 1
      end if
      lj%energy_error_tolerance = tol
    else
      lj%energy_error_tolerance = 1d-10
    end if

  end function construct_lj_type

  subroutine set_system(this, a, z, r_frac)

    use, intrinsic :: iso_fortran_env, only: error_unit
    use, non_intrinsic :: periodic_lattice_module, only: krpl => kr

    implicit none

    class(lj_type), intent(inout) :: this
    real(kr), intent(in) :: a(3, 3) ! direct lattice vectors as column vectors
    integer, intent(in) :: z(:) ! atomic numbers
    real(kr), intent(in) :: r_frac(:, :) ! fractional atomic coords.

    integer :: num_atom
    real(kr) :: rl(3)
    type(periodic_lattice_type) :: reciprocal
    real(kr) :: h(3)
    real(kr) :: h_min
    real(kr) :: r_cut
    real(kr) :: h_cut
    real(kr) :: rng(2)
    real(kr) :: delta
    real(kr) :: rij_nearest(3)
    integer :: i
    integer :: j
    integer :: m
    integer :: n

    if (size(r_frac, dim=1) /= 3 .or. size(r_frac, dim=2) /= size(z)) then
      write(error_unit, "(a, ': ', a, ': ', a)") "lj_type", "set_system", "shape of 'r_frac' is not 3 x 'size(z)'"
      stop 1
    end if

    this%lattice%vector = real(a, krpl)
    this%z = z
    this%r_frac = r_frac

    num_atom = size(z)

    if (allocated(this%a)) then
      deallocate(this%a)
    end if
    allocate(this%a(num_atom, num_atom))
    if (allocated(this%b)) then
      deallocate(this%b)
    end if
    allocate(this%b(num_atom, num_atom))
    !$omp parallel do private(m, n, i)
    do j = 1, num_atom
      n = findloc(this%param%z, z(j))
      do i = 1, num_atom
        m = findloc(this%param%z, z(i))
        this%a(i, j) = 4._kr * this%param%epsilon(m, n) * this%param%sigma(m, n) ** 12
        this%b(i, j) = -4._kr * this%param%epsilon(m, n) * this%param%sigma(m, n) ** 6
      end do
    end do
    !$omp end parallel do

    ! Eq. (33)
    rl(1) = norm(a(:, 1))
    rl(2) = norm(a(:, 2))
    rl(3) = norm(a(:, 3))
    reciprocal = this%lattice%reciprocal()
    h(1) = norm(reciprocal%vector(:, 1))
    h(2) = norm(reciprocal%vector(:, 2))
    h(3) = norm(reciprocal%vector(:, 3))
    h = h * 2 * pi
    h_min = minval(h)
    this%eta = sqrt(2._kr * minval(rl) / h_min)

    ! r_cut
    delta = huge(delta)
    !$omp parallel do reduction(min: delta), private(rij_nearest, i)
    do j = 1, num_atom
      do i = j + 1, num_atom
        rij_nearest = modulo(r_frac(:, j) - r_frac(:, i) + 0.5_kr, 1._kr) - 0.5_kr
        rij_nearest = matmul(a, rij_nearest)
        delta = min(delta, norm(rij_nearest))
      end do
    end do
    !$omp end parallel do
    if (this%energy_error_real(delta) < this%energy_error_tolerance) then
      r_cut = delta
    else
      rng = [delta, 2 * delta]
      do while (this%energy_error_tolerance < this%energy_error_real(rng(2)))
        rng = [rng(2), 2 * rng(2)]
      end do
      do while (delta < rng(2) - rng(1))
        r_cut = 0.5_kr * (rng(1) + rng(2))
        if (this%energy_error_tolerance < this%energy_error_real(r_cut)) then
          rng(1) = r_cut
        else
          rng(2) = r_cut
        end if
      end do
      r_cut = rng(2)
    end if

    ! h_cut
    delta = h_min
    if (this%energy_error_reciprocal(delta) < this%energy_error_tolerance) then
      h_cut = delta
    else
      rng = [delta, 2 * delta]
      do while (this%energy_error_tolerance < this%energy_error_reciprocal(rng(2)))
        rng = [rng(2), 2 * rng(2)]
      end do
      do while (delta < rng(2) - rng(1))
        h_cut = 0.5_kr * (rng(1) + rng(2))
        if (this%energy_error_tolerance < this%energy_error_reciprocal(h_cut)) then
          rng(1) = h_cut
        else
          rng(2) = h_cut
        end if
      end do
      h_cut = rng(2)
    end if

    if (allocated(this%cell_real)) then
      deallocate(this%cell_real)
    end if
    allocate(this%cell_real(num_atom, num_atom))

    !$omp parallel do private(i)
    do j = 1, num_atom
      do i = 1, num_atom
        this%cell_real(j, i) = this%lattice%get_cell_of_replica( &
            & real(r_frac(:, j), krpl), &
            & real(r_frac(:, i), krpl), &
            & real(r_cut, krpl))
      end do
    end do
    !$omp end parallel do

    this%cell_reciprocal = reciprocal%get_cell_of_replica( &
        & [0._kr, 0._kr, 0._kr], [0._kr, 0._kr, 0._kr], real(0.5_kr * h_cut / pi, krpl))

  end subroutine set_system

  real(kr) function energy(this)

    use, intrinsic :: iso_fortran_env, only: error_unit

    implicit none

    class(lj_type), intent(in) :: this

    integer :: num_atom
    type(periodic_lattice_type) :: reciprocal
    real(kr) :: e
    real(kr) :: rj_frac_reduced(3)
    real(kr) :: rj_frac(3)
    real(kr) :: rij(3)
    real(kr) :: dij
    real(kr) :: a
    real(kr) :: a2
    real(kr) :: h(3)
    real(kr) :: hh
    real(kr) :: b
    real(kr) :: b2
    real(kr) :: bij_cos
    integer :: i
    integer :: j
    integer :: jj
    integer :: ll

    if (.not. allocated(this%z)) then
      write(error_unit, "(a, ': ', a, ': ', a)") "lj_type", "energy", "system is not set"
      stop 1
    end if

    num_atom = size(this%z)
    reciprocal = this%lattice%reciprocal()

    energy = 0._kr

    ! r^{-12} part
    e = 0._kr
    !$omp parallel do reduction(+: e), private(rj_frac_reduced, rj_frac, rij, dij, j, jj)
    do i = 1, num_atom
      do j = 1, num_atom
        rj_frac_reduced = modulo(this%r_frac(:, j), 1._kr)
        do jj = 1, size(this%cell_real(j, i)%list, dim=2)
          if (j == i .and. all(this%cell_real(j, i)%list(:, jj) == floor(this%r_frac(:, i)))) then
            cycle
          end if
          rj_frac = this%cell_real(j, i)%list(:, jj) + rj_frac_reduced
          rij = matmul(this%lattice%vector, rj_frac - this%r_frac(:, i))
          dij = norm(rij)
          e = e + this%a(i, j) / dij ** 12
        end do
      end do
    end do
    !$omp end parallel do
    e = e * 0.5_kr
    energy = energy + e

    ! 1st term of Eq. (24), r^{-6} part
    e = 0._kr
    !$omp parallel do reduction(+: e), private(rj_frac_reduced, rj_frac, rij, dij, a, a2, j, jj)
    do i = 1, num_atom
      do j = 1, num_atom
        rj_frac_reduced = modulo(this%r_frac(:, j), 1._kr)
        do jj = 1, size(this%cell_real(j, i)%list, dim=2)
          if (j == i .and. all(this%cell_real(j, i)%list(:, jj) == floor(this%r_frac(:, i)))) then
            cycle
          end if
          rj_frac = this%cell_real(j, i)%list(:, jj) + rj_frac_reduced
          rij = matmul(this%lattice%vector, rj_frac - this%r_frac(:, i))
          dij = norm(rij)
          a = dij / this%eta
          a2 = a * a
          e = e + this%b(i, j) * (0.5_kr + (1 + 1 / a2) / a2) / a2 * exp(-a2)
        end do
      end do
    end do
    !$omp end parallel do
    e = e * 0.5_kr / this%eta ** 6
    energy = energy + e

    ! 2nd term of Eq. (24), r^{-6} part
    e = 0._kr
    !$omp parallel do reduction(+: e), private(h, hh, b, b2, rij, bij_cos, i, j)
    do ll = 1, size(this%cell_reciprocal%list, dim=2)
      if (all(this%cell_reciprocal%list(:, ll) == 0)) then
        cycle
      end if
      h = 2 * pi * matmul(reciprocal%vector, this%cell_reciprocal%list(:, ll))
      hh = norm(h)
      b = 0.5_kr * hh * this%eta
      b2 = b * b
      bij_cos = 0._kr
      do j = 1, num_atom
        do i = 1, num_atom
          rij = matmul(this%lattice%vector, this%r_frac(:, j) - this%r_frac(:, i))
          bij_cos = bij_cos + this%b(i, j) * cos(dot_product(h, rij))
        end do
      end do
      e = e + bij_cos * hh ** 3 * (sqrt_pi * erfc(b) + (0.5_kr / b2 - 1) / b * exp(-b2))
    end do
    !$omp end parallel do
    e = e * sqrt_pi ** 3 / 24._kr / this%lattice%volume()
    energy = energy + e

    ! 3rd term of Eq. (24), r^{-6} part
    energy = energy + sqrt_pi ** 3 / 6._kr / this%lattice%volume() / this%eta ** 3 * sum(this%b)

    ! 4th term of Eq. (24), r^{-6} part
    energy = energy - trace(this%b) / 12._kr / this%eta ** 6

  end function energy

  function force(this)

    use, intrinsic :: iso_fortran_env, only: error_unit

    implicit none

    class(lj_type), intent(in) :: this
    real(kr), allocatable :: force(:, :)

    integer :: num_atom
    type(periodic_lattice_type) :: reciprocal
    real(kr), allocatable :: f(:, :)
    real(kr) :: rj_frac_reduced(3)
    real(kr) :: rj_frac(3)
    real(kr) :: rij(3)
    real(kr) :: dij
    real(kr) :: a
    real(kr) :: a2
    real(kr) :: h(3)
    real(kr) :: hh
    real(kr) :: b
    real(kr) :: b2
    real(kr) :: vec(3)
    real(kr) :: bij_sin
    integer :: i
    integer :: j
    integer :: jj
    integer :: ll

    if (.not. allocated(this%z)) then
      write(error_unit, "(a, ': ', a, ': ', a)") "lj_type", "force", "system is not set"
      stop 1
    end if

    num_atom = size(this%z)
    reciprocal = this%lattice%reciprocal()

    allocate(force(3, num_atom))
    allocate(f(3, num_atom))

    force = 0._kr

    ! r^{-12} part
    f = 0._kr
    !$omp parallel do reduction(+: f), private(rj_frac_reduced, rj_frac, rij, dij, j, jj)
    do i = 1, num_atom
      do j = 1, num_atom
        rj_frac_reduced = modulo(this%r_frac(:, j), 1._kr)
        do jj = 1, size(this%cell_real(j, i)%list, dim=2)
          if (j == i .and. all(this%cell_real(j, i)%list(:, jj) == floor(this%r_frac(:, i)))) then
            cycle
          end if
          rj_frac = this%cell_real(j, i)%list(:, jj) + rj_frac_reduced
          rij = matmul(this%lattice%vector, rj_frac - this%r_frac(:, i))
          dij = norm(rij)
          f(:, i) = f(:, i) - this%a(i, j) * rij / dij ** 14
        end do
      end do
    end do
    !$omp end parallel do
    f = f * 12._kr
    force = force + f

    ! 1st term of Eq. (27), r^{-6} part
    f = 0._kr
    !$omp parallel do reduction(+: f), private(rj_frac_reduced, rj_frac, rij, dij, a, a2, j, jj)
    do i = 1, num_atom
      do j = 1, num_atom
        rj_frac_reduced = modulo(this%r_frac(:, j), 1._kr)
        do jj = 1, size(this%cell_real(j, i)%list, dim=2)
          if (j == i .and. all(this%cell_real(j, i)%list(:, jj) == floor(this%r_frac(:, i)))) then
            cycle
          end if
          rj_frac = this%cell_real(j, i)%list(:, jj) + rj_frac_reduced
          rij = matmul(this%lattice%vector, rj_frac - this%r_frac(:, i))
          dij = norm(rij)
          a = dij / this%eta
          a2 = a * a
          f(:, i) = f(:, i) - this%b(i, j) * (1 + 3 * (1 + 2 * (1 + 1 / a2) / a2) / a2) / a2 * exp(-a2) * rij
        end do
      end do
    end do
    !$omp end parallel do
    f = f / this%eta ** 8
    force = force + f

    ! 2nd term of Eq. (27), r^{-6} part
    f = 0._kr
    !$omp parallel do reduction(+: f), private(h, hh, b, b2, vec, rij, bij_sin, i, j)
    do ll = 1, size(this%cell_reciprocal%list, dim=2)
      if (all(this%cell_reciprocal%list(:, ll) == 0)) then
        cycle
      end if
      h = 2 * pi * matmul(reciprocal%vector, this%cell_reciprocal%list(:, ll))
      hh = norm(h)
      b = 0.5_kr * hh * this%eta
      b2 = b * b
      vec = hh ** 3 * (sqrt_pi * erfc(b) + (0.5_kr / b2 - 1) / b * exp(-b2)) * h
      do i = 1, num_atom
        bij_sin = 0._kr
        do j = 1, num_atom
          rij = matmul(this%lattice%vector, this%r_frac(:, j) - this%r_frac(:, i))
          bij_sin = bij_sin - this%b(i, j) * sin(dot_product(h, rij))
        end do
        f(:, i) = f(:, i) + bij_sin * vec
      end do
    end do
    !$omp end parallel do
    f = f * sqrt_pi ** 3 / 12._kr / this%lattice%volume()
    force = force + f

    deallocate(f)

  end function force

  function stress(this)

    use, intrinsic :: iso_fortran_env, only: error_unit

    implicit none

    class(lj_type), intent(in) :: this
    real(kr) :: stress(3, 3)

    integer :: num_atom
    type(periodic_lattice_type) :: reciprocal
    real(kr) :: s(3, 3)
    real(kr) :: rj_frac_reduced(3)
    real(kr) :: rj_frac(3)
    real(kr) :: rij(3)
    real(kr) :: dij
    real(kr) :: a
    real(kr) :: a2
    real(kr) :: h(3)
    real(kr) :: hh
    real(kr) :: b
    real(kr) :: b2
    real(kr) :: c
    real(kr) :: mat(3, 3)
    real(kr) :: bij_cos
    integer :: i
    integer :: j
    integer :: jj
    integer :: ll

    if (.not. allocated(this%z)) then
      write(error_unit, "(a, ': ', a, ': ', a)") "lj_type", "stress", "system is not set"
      stop 1
    end if

    num_atom = size(this%z)
    reciprocal = this%lattice%reciprocal()

    stress = 0._kr

    ! r^{-12} part
    s = 0._kr
    !$omp parallel do reduction(+: s), private(rj_frac_reduced, rj_frac, rij, dij, j, jj)
    do i = 1, num_atom
      do j = 1, num_atom
        rj_frac_reduced = modulo(this%r_frac(:, j), 1._kr)
        do jj = 1, size(this%cell_real(j, i)%list, dim=2)
          if (j == i .and. all(this%cell_real(j, i)%list(:, jj) == floor(this%r_frac(:, i)))) then
            cycle
          end if
          rj_frac = this%cell_real(j, i)%list(:, jj) + rj_frac_reduced
          rij = matmul(this%lattice%vector, rj_frac - this%r_frac(:, i))
          dij = norm(rij)
          s = s + this%a(i, j) / dij ** 14 * outer_product(rij, rij)
        end do
      end do
    end do
    !$omp end parallel do
    s = s * 6._kr
    stress = stress + s

    ! 1st term of Eq. (28), r^{-6} part
    s = 0._kr
    !$omp parallel do reduction(+: s), private(rj_frac_reduced, rj_frac, rij, dij, a, a2, j, jj)
    do i = 1, num_atom
      do j = 1, num_atom
        rj_frac_reduced = modulo(this%r_frac(:, j), 1._kr)
        do jj = 1, size(this%cell_real(j, i)%list, dim=2)
          if (j == i .and. all(this%cell_real(j, i)%list(:, jj) == floor(this%r_frac(:, i)))) then
            cycle
          end if
          rj_frac = this%cell_real(j, i)%list(:, jj) + rj_frac_reduced
          rij = matmul(this%lattice%vector, rj_frac - this%r_frac(:, i))
          dij = norm(rij)
          a = dij / this%eta
          a2 = a * a
          s = s + this%b(i, j) * (1 + 3 * (1 + 2 * (1 + 1 / a2) / a2) / a2) / a2 * exp(-a2) * outer_product(rij, rij)
        end do
      end do
    end do
    !$omp end parallel do
    s = s * 0.5_kr / this%eta ** 8
    stress = stress + s

    ! 2nd term of Eq. (28), r^{-6} part
    s = 0._kr
    !$omp parallel do reduction(+: s), private(h, hh, b, b2, c, mat, rij, bij_cos, i, j)
    do ll = 1, size(this%cell_reciprocal%list, dim=2)
      if (all(this%cell_reciprocal%list(:, ll) == 0)) then
        cycle
      end if
      h = 2 * pi * matmul(reciprocal%vector, this%cell_reciprocal%list(:, ll))
      hh = norm(h)
      b = 0.5_kr * hh * this%eta
      b2 = b * b
      c = sqrt_pi * erfc(b) - exp(-b2) / b
      mat = 0._kr
      mat(1, 1) = hh ** 3 * (c + 0.5_kr * exp(-b2) / b ** 3)
      mat(2, 2) = mat(1, 1)
      mat(3, 3) = mat(1, 1)
      mat = mat + 3._kr * hh * c * outer_product(h, h)
      bij_cos = 0._kr
      do j = 1, num_atom
        do i = 1, num_atom
          rij = matmul(this%lattice%vector, this%r_frac(:, j) - this%r_frac(:, i))
          bij_cos = bij_cos + this%b(i, j) * sin(dot_product(h, rij))
        end do
      end do
      s = s + bij_cos * mat
    end do
    !$omp end parallel do
    s = s * 0.5_kr * sqrt_pi ** 3 / 12.0_kr / this%lattice%volume()
    stress = stress + s

    ! 3rd term of Eq. (28), r^{-6} part
    s(1, 1) = 0.5_kr * sqrt_pi ** 3 / 3._kr / this%eta ** 3 / this%lattice%volume() * sum(this%b)
    stress(1, 1) = stress(1, 1) + s(1, 1)
    stress(2, 2) = stress(2, 2) + s(1, 1)
    stress(3, 3) = stress(3, 3) + s(1, 1)

    stress = stress / this%lattice%volume()

  end function stress

  pure real(kr) function energy_error_real(this, r_cut)

    implicit none

    class(lj_type), intent(in) :: this
    real(kr), intent(in) :: r_cut

    integer :: num_atom
    real(kr) :: n2a
    real(kr) :: n2b
    real(kr) :: r_cut_over_eta
    integer :: i
    integer :: j

    num_atom = size(this%z)

    ! Eq. (45)
    n2a = sum(this%a)
    energy_error_real = 2 * pi * n2a / 9._kr / this%lattice%volume() / r_cut ** 9

    ! Eq. (40)
    n2b = sum(abs(this%b))
    r_cut_over_eta = r_cut / this%eta
    energy_error_real = energy_error_real - &
        & sqrt_pi ** 3 * n2b / this%lattice%volume() * this%eta &
        & * (1 + (1 + 0.5_kr * r_cut_over_eta * r_cut_over_eta) * r_cut_over_eta * r_cut_over_eta) / r_cut ** 4  &
        & * erfc(r_cut_over_eta)

    energy_error_real = abs(energy_error_real)

  end function energy_error_real

  ! Eq. (42)
  pure real(kr) function energy_error_reciprocal(this, h_cut)

    implicit none

    class(lj_type), intent(in) :: this
    real(kr), intent(in) :: h_cut

    real(kr) :: n2b
    real(kr) :: h_cut_eta

    n2b = sum(abs(this%b))
    h_cut_eta = h_cut * this%eta

    energy_error_reciprocal = n2b / 6._kr / sqrt_pi / this%eta ** 6 &
        & * (h_cut_eta * exp(-0.25_kr * h_cut_eta * h_cut_eta) + sqrt_pi * erfc(0.5_kr * h_cut_eta))

  end function energy_error_reciprocal

  type(lj_parameter_type) function make_parameter_lorentz_berthelot(z, epsilon, sigma) result(param)

    use, intrinsic :: iso_fortran_env, only: error_unit

    implicit none

    integer, intent(in) :: z(:)
    real(kr), intent(in) :: epsilon(:)
    real(kr), intent(in) :: sigma(:)

    if (size(epsilon) /= size(z) .or. size(sigma) /= size(z)) then
      write(error_unit, "(a, ': ', a)") "make_parameter_lorentz_berthelot", "size mismatch for 'z', 'epsilon', and 'sigma'"
      stop 1
    end if

    param%z = z
    param%epsilon = sqrt(outer_product(epsilon, epsilon))
    param%sigma = 0.5_kr * outer_sum(sigma, sigma)

  end function make_parameter_lorentz_berthelot

  type(lj_parameter_type) function make_parameter_kong(z, epsilon, sigma) result(param)

    use, intrinsic :: iso_fortran_env, only: error_unit

    implicit none

    integer, intent(in) :: z(:)
    real(kr), intent(in) :: epsilon(:)
    real(kr), intent(in) :: sigma(:)

    real(kr), allocatable :: vec(:)
    real(kr), allocatable :: mat1(:, :)
    real(kr), allocatable :: mat2(:, :)

    if (size(epsilon) /= size(z) .or. size(sigma) /= size(z)) then
      write(error_unit, "(a, ': ', a)") "make_parameter_kong", "size mismatch for 'z', 'epsilon', and 'sigma'"
      stop 1
    end if

    vec = (epsilon * sigma ** 12) ** (1._kr / 13._kr)
    mat1 = (0.5_kr * outer_sum(vec, vec)) ** 13

    vec = epsilon * sigma ** 6
    mat2 = outer_product(vec, vec)

    param%z = z
    param%epsilon = mat2 / mat1
    param%sigma = (mat1 / sqrt(mat2)) ** (1._kr / 6._kr)

  end function make_parameter_kong

  type(lj_parameter_type) function make_parameter_waldman_hagler(z, epsilon, sigma) result(param)

    use, intrinsic :: iso_fortran_env, only: error_unit

    implicit none

    integer, intent(in) :: z(:)
    real(kr), intent(in) :: epsilon(:)
    real(kr), intent(in) :: sigma(:)

    real(kr), allocatable :: vec(:)
    real(kr), allocatable :: mat(:, :)

    if (size(epsilon) /= size(z) .or. size(sigma) /= size(z)) then
      write(error_unit, "(a, ': ', a)") "make_parameter_kong", "size mismatch for 'z', 'epsilon', and 'sigma'"
      stop 1
    end if

    vec = sigma ** 3
    mat = 0.5_kr * outer_sum(vec * vec, vec * vec)

    param%z = z
    param%epsilon = outer_product(vec, vec) / mat * sqrt(outer_product(epsilon, epsilon))
    param%sigma = mat ** (1._kr / 6._kr)

  end function make_parameter_waldman_hagler

  pure real(kr) function norm(v)

    implicit none

    real(kr), intent(in) :: v(:)

    norm = sqrt(sum(v * v))

  end function norm

  pure real(kr) function trace(m)

    implicit none

    real(kr), intent(in) :: m(:, :)

    integer :: i

    trace = 0._kr
    do i = 1, minval(shape(m))
      trace = trace + m(i, i)
    end do

  end function trace

  pure function outer_sum(a, b)

    implicit none

    real(kr), intent(in) :: a(:)
    real(kr), intent(in) :: b(:)
    real(kr), allocatable :: outer_sum(:, :)

    integer :: j

    allocate(outer_sum(size(a), size(b)))
    do j = 1, size(b)
      outer_sum(:, j) = a + b(j)
    end do

  end function outer_sum

  pure function outer_product(a, b)

    implicit none

    real(kr), intent(in) :: a(:)
    real(kr), intent(in) :: b(:)
    real(kr), allocatable :: outer_product(:, :)

    integer :: j

    allocate(outer_product(size(a), size(b)))
    do j = 1, size(b)
      outer_product(:, j) = a * b(j)
    end do

  end function outer_product

  pure integer function findloc(array, value)

    implicit none

    integer, intent(in) :: array(:)
    integer, intent(in) :: value

    integer :: i

    findloc = 0
    do i = 1, size(array)
      if (array(i) == value) then
        findloc = i
        exit
      end if
    end do

  end function findloc

end module lj_module
