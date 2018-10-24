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

module m_jmatgen

  implicit none

  private

  integer, parameter, public :: kr = kind(0d0)

  public :: initJmatgen, calcJmatgen

contains

  subroutine initJmatgen(config_energy, config_force)

    character(len=*), intent(in) :: config_energy
    character(len=*), intent(in) :: config_force

    ! do nothing

  end subroutine initJmatgen

  subroutine calcJmatgen(natom, cell, atom, coord, energy, force, virial)

    use, intrinsic :: iso_fortran_env, only: error_unit

    integer, intent(in) :: natom
    real(kr), intent(in) :: cell(3, 3)
    character(len=*), intent(in) :: atom(natom)
    real(kr), intent(in) :: coord(3, natom)
    real(kr), intent(out) :: energy
    real(kr), intent(inout) :: force(3, natom)
    real(kr), intent(inout) :: virial(3, 3)

    write(error_unit, "(a)") "not compiled to use Jmatgen"
    stop 1

  end subroutine calcJmatgen

end module m_jmatgen
