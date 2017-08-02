module periodic_lattice_module

  implicit none

  private

  integer, parameter, public :: kr = kind(0d0)

  type, public :: periodic_lattice_type

    private

    real(kr), public :: direct_lattice(3, 3) ! direct lattice vectors as column vectors

  contains

    private

    procedure, public :: volume
    procedure, public :: reciprocal_lattice
    procedure, public :: get_cell_of_replica

  end type periodic_lattice_type

contains

  pure real(kr) function volume(this)

    implicit none

    class(periodic_lattice_type), intent(in) :: this

    volume = dot_product(cross_product(this%direct_lattice(:, 1), this%direct_lattice(:, 2)), this%direct_lattice(:, 3))

  end function volume

  pure function reciprocal_lattice(this)

    implicit none

    class(periodic_lattice_type), intent(in) :: this
    real(kr) :: reciprocal_lattice(3, 3)

    ! without 2 pi
    reciprocal_lattice(:, 1) = cross_product(this%direct_lattice(:, 2), this%direct_lattice(:, 3))
    reciprocal_lattice(:, 2) = cross_product(this%direct_lattice(:, 3), this%direct_lattice(:, 1))
    reciprocal_lattice(:, 3) = cross_product(this%direct_lattice(:, 1), this%direct_lattice(:, 2))
    reciprocal_lattice = reciprocal_lattice / this%volume()

  end function reciprocal_lattice

  function get_cell_of_replica(this, r_frac, origin_frac, cutoff) result(cell_replica)

    use, intrinsic :: iso_fortran_env, only: error_unit

    implicit none

    class(periodic_lattice_type), intent(in) :: this
    real(kr), intent(in) :: r_frac(3) ! target point in fractional coords.
    real(kr), intent(in) :: origin_frac(3) ! origin of a searched region in fractional coords.
    real(kr), intent(in) :: cutoff ! cutoff of a searched region
    integer, allocatable :: cell_replica(:, :)

    real(kr) :: b(3, 3)
    real(kr) :: cutoff_cell(3)
    integer :: cell_min(3)
    integer :: cell_max(3)
    integer :: cell(3)
    real(kr) :: r_frac_reduced(3)
    real(kr) :: r_replica_frac(3)
    real(kr) :: r_replica(3)
    real(kr) :: origin(3)
    real(kr) :: d
    integer :: num_replica
    integer :: n1
    integer :: n2
    integer :: n3

    ! 1/8 size of a supercell circumscribing a sphere of radius 'r_cut'
    b = this%reciprocal_lattice()
    cutoff_cell(1) = norm(b(:, 1)) * cutoff
    cutoff_cell(2) = norm(b(:, 2)) * cutoff
    cutoff_cell(3) = norm(b(:, 3)) * cutoff

    if (any(exceed_integer(origin_frac - cutoff_cell)) &
        & .or. any(exceed_integer(origin_frac + cutoff_cell))) then
      write(error_unit, "(a, ': ', a, ': ', a)") "periodic_lattice_type", "get_cell_of_replica", "exceed the range of 'integer'"
      stop 1
    end if
    cell_min = floor(origin_frac - cutoff_cell)
    cell_max = floor(origin_frac + cutoff_cell)

    r_frac_reduced = modulo(r_frac, 1._kr)
    origin = matmul(this%direct_lattice, origin_frac)

    allocate(cell_replica(3, product(cell_max - cell_min + 1)))
    num_replica = 0

    do n3 = cell_min(3), cell_max(3)
      cell(3) = n3
      do n2 = cell_min(2), cell_max(2)
        cell(2) = n2
        do n1 = cell_min(1), cell_max(1)
          cell(1) = n1

          r_replica_frac = cell + r_frac_reduced
          r_replica = matmul(this%direct_lattice, r_replica_frac)
          d = norm(r_replica - origin)
          if (d > cutoff) then
            cycle
          end if

          cell_replica(:, num_replica + 1) = cell
          num_replica = num_replica + 1

        end do
      end do
    end do

    cell_replica = cell_replica(:, 1:num_replica)

  end function get_cell_of_replica

  elemental logical function exceed_integer(r)

    implicit none

    real(kr), intent(in) :: r

    integer, parameter :: integer_max = huge(0)

    exceed_integer = r < -real(integer_max, kr) - 1._kr .or. real(integer_max, kr) < r

  end function exceed_integer

  pure real(kr) function norm(v)

    real(kr), intent(in) :: v(:)

    norm = sqrt(sum(v ** 2))

  end function norm

  pure function cross_product(v1, v2)

    implicit none

    real(kr), intent(in) :: v1(3)
    real(kr), intent(in) :: v2(3)
    real(kr) :: cross_product(3)

    cross_product(1) = v1(2) * v2(3) - v1(3) * v2(2)
    cross_product(2) = v1(3) * v2(1) - v1(1) * v2(3)
    cross_product(3) = v1(1) * v2(2) - v1(2) * v2(1)

  end function cross_product

end module periodic_lattice_module
