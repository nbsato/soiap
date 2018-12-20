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

subroutine show_argument()
  implicit none
  write(6,*)
  write(6,*)'argument:'
  write(6,*)'thisprogram inputfilename'
  write(6,*)
end subroutine show_argument

subroutine input
  use keyvalue
  use paramlist
  use cif_module
  use lj_module, only: &
      & lj_parameter_type, &
      & make_parameter_lorentz_berthelot, &
      & make_parameter_kong, &
      & make_parameter_waldman_hagler
  use m_jmatgen, only: initJmatgen

  implicit none

  integer :: ifin,ret,ilength,icolumn,iatompos,itime,mscale,isymmetry
  integer :: i,j,itmp
  real(8) :: uvin(3,3),vtmp(3),aw,volscale,rscale,rtmp
  character :: ctmp*2

  character(120):: inputfilename
  character(120):: ciffilename
  character(256):: config_jmatgen_energy
  character(256):: config_jmatgen_force
  integer:: narg

  narg=command_argument_count()
  if (narg/=1) then 
     call show_argument()
     stop 100
  endif
  call get_command_argument(1,inputfilename)
  write(6,*)'input=',trim(inputfilename)

  QMD%is_symmetrized=.false.
  call getkeyvalue(inputfilename,"crystal",ciffilename,default="",unit=ifin,status=ret)
  if( ciffilename /= "" ) then
     call CIF_load(ciffilename)
     call getkeyvalue(inputfilename,"symmetry",isymmetry,default=0)
     QMD%is_symmetrized=isymmetry==1
  else
     ! unit cell
     call getkeyvalue(inputfilename,"unit_vec",unit=ifin,status=ret)
     write(*,*)'unit_vec (input)'
     do i=1,3
        read(ifin,*)ctmp,(uvin(i,j),j=1,3)
        write(*,"(3f23.16)")(uvin(i,j),j=1,3)
     enddo
     close(ifin)

     call getkeyvalue(inputfilename,"unit_vec_column",icolumn,default=1)
     if (icolumn==2) then
        QMD%uv=transpose(uvin)
        write(*,*)'unit vector (column/raw): raw'
     else
        QMD%uv=uvin
        write(*,*)'unit vector (column/raw): column'
     endif

     call getkeyvalue(inputfilename,"unit_length",ilength,default=1)
     if (ilength==1) QMD%uv=QMD%uv/bohr

     call getkeyvalue(inputfilename,"volume_scale",volscale,default=1.d0)
     rscale=volscale**(1d0/3d0)
     if (abs(rscale-1.d0).gt.1.d-6) then 
        write(*,*)'volume_scale =',volscale
        QMD%uv=QMD%uv*rscale
     endif

     ! u1=QMD%uv(1:3,1)
     ! u2=QMD%uv(1:3,2)
     ! u3=QMD%uv(1:3,3)

     call cross_x(QMD%uv(:,1),QMD%uv(:,2),QMD%bv(:,3))
     call cross_x(QMD%uv(:,2),QMD%uv(:,3),QMD%bv(:,1))
     call cross_x(QMD%uv(:,3),QMD%uv(:,1),QMD%bv(:,2))
     QMD%omega=dot_product(QMD%bv(:,3),QMD%uv(:,3))
     QMD%omegai=1.d0/QMD%omega
     QMD%bv=QMD%bv*QMD%omegai
     ! dot_product(ai,bj) = delta(i,j) <-- not 2pi*delta(i,j)

     ! informtion of atoms
     call getkeyvalue(inputfilename,"number_atom",QMD%natom,default=0)
     allocate(QMD%ra(3,QMD%natom),QMD%rr(3,QMD%natom),QMD%iposfix(QMD%natom))
     allocate(QMD%frc(3,QMD%natom),QMD%vrr(3,QMD%natom))
     QMD%frc=0.d0
     QMD%vrr=0.d0

     call getkeyvalue(inputfilename,"number_element",QMD%nkatm,default=0)
     allocate(QMD%zatm(QMD%natom),QMD%mass(QMD%natom),QMD%mfac(QMD%natom))
     !  QMD%mconv=(1.6605655d-27/9.109534d-31)*0.5d0 ! Rydberg atomic units
     QMD%mconv=(1.6605655d-27/9.109534d-31) ! Hartree atomic units

     call getkeyvalue(inputfilename,"atom_pos",iatompos,default=1)

     call getkeyvalue(inputfilename,"atom_list",unit=ifin,status=ret)
     do i=1,QMD%natom
        read(ifin,*)QMD%zatm(i),vtmp(1:3),QMD%iposfix(i)
        if (iatompos==1) then
           QMD%rr(:,i)=vtmp
           QMD%ra(:,i)=matmul(QMD%uv,vtmp)
        else
           if (ilength==1) vtmp=vtmp/bohr
           if (abs(rscale-1.d0).gt.1.d-6) then 
              !           write(*,*)'volume_scale =',volscale
              vtmp=vtmp*rscale
           endif

           QMD%ra(:,i)=vtmp
           QMD%rr(:,i)=matmul(vtmp,QMD%bv)
        endif
     enddo
     close(ifin)
  end if ! ciffilename

  write(*,*)'unit vectors [Bohr]'
  write(*,"('u1 =',3f23.16)")QMD%uv(1:3,1)
  write(*,"('u2 =',3f23.16)")QMD%uv(1:3,2)
  write(*,"('u3 =',3f23.16)")QMD%uv(1:3,3)
  write(*,*)'QMD%omega =',QMD%omega

  write(*,*)'number of atoms =',QMD%natom
  write(*,*)'number of elements =',QMD%nkatm

  write(*,*)'atom_list: lattice coordinate'
  do i=1,QMD%natom
     write(*,"(i4,3f23.16,i4)")QMD%zatm(i),(QMD%rr(j,i),j=1,3), &
          QMD%iposfix(i)
  enddo

  write(*,*)'atom_list: Cartesian coordinate [Bohr]'
  do i=1,QMD%natom
     write(*,"(i4,3f23.16,i4)")QMD%zatm(i),(QMD%ra(j,i),j=1,3), &
          QMD%iposfix(i)
  enddo

  if (QMD%is_symmetrized) write(*,*)'symmetry on'

! optimization
  call getkeyvalue(inputfilename,"md_mode",QMD%imd,default=4)
  write(*,*)'md_mode =',QMD%imd

  call getkeyvalue(inputfilename,"number_max_relax",QMD%nloopa,default=10)
  write(*,*)'number_max_relax =',QMD%nloopa

  call getkeyvalue(inputfilename,"md_mode_cell",QMD%imdc,default=2)
  write(*,*)'md_mode_cell =',QMD%imdc

  call getkeyvalue(inputfilename,"number_max_relax_cell",QMD%nloopc,default=0)
  write(*,*)'number_max_relax_cell =',QMD%nloopc

  call getkeyvalue(inputfilename,"external_stress_v",QMD%extstrs,3,default=(/0.0d0,0.0d0,0.0d0/))
  write(*,"(a17,3f12.6)")'external_stress_v [GPa] =',QMD%extstrs
  QMD%extstrs=QMD%extstrs/29421.010901602753 ! GPa to Ha/Bohr^3

  call getkeyvalue(inputfilename,"time_step",QMD%tstep,default=1.d0)
  write(*,*)'time_step (input) =',QMD%tstep

  call getkeyvalue(inputfilename,"unit_time",itime,default=1)
  if (itime==1) QMD%tstep=QMD%tstep*100/2.418884326505
  write(*,*)'time_step [a.u.] =',QMD%tstep
  QMD%tstep0=QMD%tstep
  QMD%tstepc=QMD%tstep

  call getkeyvalue(inputfilename,"th_force",QMD%fth,default=5.d-5)
  write(*,*)'th_force =',QMD%fth

  call getkeyvalue(inputfilename,"th_stress",QMD%sth,default=5.d-7)
  write(*,*)'th_stress =',QMD%sth

  if (QMD%imdc==3) then ! RFC5
     call getkeyvalue(inputfilename,"max_displacement",QMD%rdmax,default=0.10d0)
  else
     call getkeyvalue(inputfilename,"max_displacement",QMD%rdmax,default=0.25d0)
  endif
  write(*,*)'max_displacement [Bohr] =',QMD%rdmax

  call getkeyvalue(inputfilename,"mass_scale",mscale,default=0)
  if (mscale==1) write(*,*)'mass_scale on'

  call getkeyvalue(inputfilename,"mass_cell",QMD%mcell,default=5.d-4)
  write(*,*)'mass_cell =',QMD%mcell

! atom QMD%mass
  if (QMD%imd<=0.or.mscale==1) then
     do i=1,QMD%natom
        QMD%mass(i)=12.d0*QMD%mconv
        QMD%mfac(i)=1.d0
     enddo   
  else   
     do i=1,QMD%natom
        call massset(QMD%zatm(i),aw)
        QMD%mass(i)=aw*QMD%mconv
        QMD%mfac(i)=1.d0
     enddo   
  endif

! force field
  call getkeyvalue(inputfilename,"force_field",QMD%ifrcf,default=0)
  write(*,*)'force_field =',QMD%ifrcf

  select case (QMD%ifrcf)
    case (5) ! Jmatgen
      call getkeyvalue(inputfilename,"config_jmatgen_energy",config_jmatgen_energy)
      call getkeyvalue(inputfilename,"config_jmatgen_force",config_jmatgen_force,default="")
      write(*,*)'config_jmatgen_energy=',trim(config_jmatgen_energy)
      write(*,*)'config_jmatgen_force=',trim(config_jmatgen_force)
      call initJmatgen(config_jmatgen_energy,config_jmatgen_force)
    case (6) ! Lennard-Jones potential
      allocate(QMD%lj_parameter%z(QMD%nkatm))
      call getkeyvalue(inputfilename,"lj_z",QMD%lj_parameter%z(:),size(QMD%lj_parameter%z))
      write(*,*)'lj_z=',QMD%lj_parameter%z
      allocate(QMD%lj_parameter%epsilon(QMD%nkatm, 1))
      call getkeyvalue(inputfilename,"lj_epsilon",QMD%lj_parameter%epsilon(:, 1),size(QMD%lj_parameter%epsilon))
      write(*,*)'lj_epsilon=',QMD%lj_parameter%epsilon
      allocate(QMD%lj_parameter%sigma(QMD%nkatm, 1))
      call getkeyvalue(inputfilename,"lj_sigma",QMD%lj_parameter%sigma(:, 1),size(QMD%lj_parameter%sigma))
      write(*,*)'lj_sigma=',QMD%lj_parameter%sigma
      call getkeyvalue(inputfilename,"lj_combination_rules",itmp)
      write(*,*)'lj_combination_rules=',itmp
      select case (itmp)
        case (1)
          QMD%lj_parameter=make_parameter_lorentz_berthelot(QMD%lj_parameter%z,QMD%lj_parameter%epsilon(:,1),QMD%lj_parameter%sigma(:,1))
        case (2)
          QMD%lj_parameter=make_parameter_kong(QMD%lj_parameter%z,QMD%lj_parameter%epsilon(:,1),QMD%lj_parameter%sigma(:,1))
        case (3)
          QMD%lj_parameter=make_parameter_waldman_hagler(QMD%lj_parameter%z,QMD%lj_parameter%epsilon(:,1),QMD%lj_parameter%sigma(:,1))
        case default
          stop 'invalid lj_combination_rule'
      end select
  end select

end subroutine input      

! *********************************************************************
!  2005.04.27  T.Tamura
!      subroutine massset(zatom,aw)
      subroutine massset(izatom,aw)
! *********************************************************************
      implicit none
!
      integer izatom
      real*8 zatom,aw
!
! commented out, 150610, TM
!      izatom=idint(zatom)
      aw=0.0d0
!
      if(izatom.eq.1) aw=1.0079d0
      if(izatom.eq.2) aw=4.003d0
      if(izatom.eq.3) aw=6.941d0
      if(izatom.eq.4) aw=9.012d0
      if(izatom.eq.5) aw=10.811d0
      if(izatom.eq.6) aw=12.011d0
      if(izatom.eq.7) aw=14.007d0
      if(izatom.eq.8) aw=15.999d0
      if(izatom.eq.9) aw=18.998d0
      if(izatom.eq.10) aw=20.180d0
      if(izatom.eq.11) aw=11.990d0
      if(izatom.eq.12) aw=24.305d0
      if(izatom.eq.13) aw=26.982d0
      if(izatom.eq.14) aw=28.086d0
      if(izatom.eq.15) aw=30.974d0
      if(izatom.eq.16) aw=32.066d0
      if(izatom.eq.17) aw=35.453d0
      if(izatom.eq.18) aw=39.948d0
      if(izatom.eq.19) aw=39.098d0
      if(izatom.eq.20) aw=40.078d0
      if(izatom.eq.21) aw=44.956d0
      if(izatom.eq.22) aw=47.88d0
      if(izatom.eq.23) aw=50.942d0
      if(izatom.eq.24) aw=51.996d0
      if(izatom.eq.25) aw=54.938d0
      if(izatom.eq.26) aw=55.847d0
      if(izatom.eq.27) aw=58.933d0
      if(izatom.eq.28) aw=58.69d0
      if(izatom.eq.29) aw=63.546d0
      if(izatom.eq.30) aw=65.39d0
      if(izatom.eq.31) aw=69.723d0
      if(izatom.eq.32) aw=72.61d0
      if(izatom.eq.33) aw=74.922d0
      if(izatom.eq.34) aw=78.96d0
      if(izatom.eq.35) aw=79.904d0
      if(izatom.eq.36) aw=83.80d0
      if(izatom.eq.37) aw=85.468d0
      if(izatom.eq.38) aw=87.62d0
      if(izatom.eq.39) aw=88.906d0
      if(izatom.eq.40) aw=91.224d0
      if(izatom.eq.41) aw=92.906d0
      if(izatom.eq.42) aw=95.94d0
      if(izatom.eq.43) aw=98.0d0
      if(izatom.eq.44) aw=101.07d0
      if(izatom.eq.45) aw=102.91d0
      if(izatom.eq.46) aw=106.42d0
      if(izatom.eq.47) aw=107.87d0
      if(izatom.eq.48) aw=112.41d0
      if(izatom.eq.49) aw=114.82d0
      if(izatom.eq.50) aw=118.71d0
      if(izatom.eq.51) aw=121.75d0
      if(izatom.eq.52) aw=127.6d0
      if(izatom.eq.53) aw=126.9d0
      if(izatom.eq.54) aw=131.29d0
      if(izatom.eq.55) aw=132.91d0
      if(izatom.eq.56) aw=137.33d0
      if(izatom.eq.57) aw=138.91d0
      if(izatom.eq.58) aw=140.12d0
      if(izatom.eq.59) aw=140.91d0
      if(izatom.eq.60) aw=144.24d0
      if(izatom.eq.61) aw=145.0d0
      if(izatom.eq.62) aw=150.36d0
      if(izatom.eq.63) aw=151.97d0
      if(izatom.eq.64) aw=157.25d0
      if(izatom.eq.65) aw=158.93d0
      if(izatom.eq.66) aw=162.5d0
      if(izatom.eq.67) aw=164.93d0
      if(izatom.eq.68) aw=167.26d0
      if(izatom.eq.69) aw=168.93d0
      if(izatom.eq.70) aw=173.04d0
      if(izatom.eq.71) aw=174.97d0
      if(izatom.eq.72) aw=178.49d0
      if(izatom.eq.73) aw=180.95d0
      if(izatom.eq.74) aw=183.85d0
      if(izatom.eq.75) aw=186.21d0
      if(izatom.eq.76) aw=190.2d0
      if(izatom.eq.77) aw=192.22d0
      if(izatom.eq.78) aw=195.08d0
      if(izatom.eq.79) aw=196.97d0
      if(izatom.eq.80) aw=200.59d0
      if(izatom.eq.81) aw=204.38d0
      if(izatom.eq.82) aw=207.2d0
      if(izatom.eq.83) aw=208.98d0
      if(izatom.eq.84) aw=209.0d0
      if(izatom.eq.85) aw=210.0d0
      if(izatom.eq.86) aw=222.0d0
      if(izatom.eq.87) aw=223.0d0
      if(izatom.eq.88) aw=226.0d0
      if(izatom.eq.89) aw=227.0d0
      if(izatom.eq.90) aw=232.04d0
      if(izatom.eq.91) aw=231.04d0
      if(izatom.eq.92) aw=238.03d0
      if(izatom.eq.93) aw=237.0d0
      if(izatom.eq.94) aw=244.0d0
      if(izatom.eq.95) aw=243.0d0
      if(izatom.eq.96) aw=247.0d0
      if(izatom.eq.97) aw=247.0d0
      if(izatom.eq.98) aw=251.0d0
      if(izatom.eq.99) aw=252.0d0
      if(izatom.eq.100) aw=257.0d0
      if(izatom.eq.101) aw=258.0d0
      if(izatom.eq.102) aw=259.0d0
      if(izatom.eq.103) aw=260.0d0
!
      if(aw.ne.0.0d0) go to 100
      stop 'no data for aw'
!
  100 continue
!
      return

      end
!--------1---------2---------3---------4---------5---------6---------7--
