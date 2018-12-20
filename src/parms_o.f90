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

module paramlist

  use lj_module, only: lj_parameter_type

  real*8 :: bohr,kcalpermol,hartree,tol
  parameter (bohr=0.529177249d0, kcalpermol=627.503d0, hartree=27.2116d0)
  parameter (tol=1.d-12)
type t_qmd

  integer :: natom,nloopa,nloopc,loopa,loopc,imd,imdc,isd_cell
  integer :: nkatm
  integer :: ifrcf
  real*8 :: uv(3,3),uvo(3,3,2),vuv(3,3),bv(3,3),omega,omegai,mconv
  real*8 :: tote
  real*8 :: strs(3,3),strso(3,3,2),extstrs(3),fth,fmax,sth,smax,rdmax
  real*8 :: tstep,tstep0,mcell,guv(3,3,0:2),duv0(3,3)
  logical :: is_symmetrized
  integer,allocatable :: iposfix(:),zatm(:)
  real*8,allocatable :: ra(:,:),rr(:,:),mass(:),mfac(:)
  real*8,allocatable :: frc(:,:),vrr(:,:)

  type(lj_parameter_type) :: lj_parameter

! lattice_fire
  integer :: fire_nminc ! for FIRE
  real*8 :: fire_fincc,fire_fdecc,fire_alpc,fire_alp0c,fire_falpc,fire_dtmaxc
  real*8 :: tstepc

! atomrelax_o.f90
  integer :: npstv,fire_nmin ! for FIRE
  real*8 :: fire_finc,fire_fdec,fire_alp,fire_alp0,fire_falp,fire_dtmax

end type
type(t_qmd):: QMD
end module paramlist
