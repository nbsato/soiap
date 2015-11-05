program main
      use m_showtime
  use paramlist
  implicit none

  integer i,j,ia,ika
      type(t_showtime)::  timer

  open(901,file='structure.out',form='formatted')

  call input
  uvo(:,:,1)=uv
  strso(:,:,1)=strs
  vuv=0.d0

! debug:
  if (imd/=0) stop 'debug: imd should be 0'
  if (imd/=0.and.imdc/=2) stop 'debug: imdc should be 0 or 2'

  do loopc=1,nloopc ! relax unit vectors
     if (loopc>1) then 
        call updt_cell
     endif   
     call output_struc(901)

     do loopa=1,nloopa ! relax internal coordinates
        call updt_coord
      call timer%start()
        call tote_frc_strs
      call timer%stop()
      call timer%show("tote_frc_strs:")
        if (fmax<fth) then 
           write(*,'("frc converged. loopc, loopa",2i5,e12.4)')loopc,loopa,fmax
           exit
        endif   
     enddo ! loopa

     do i=1,3
        strs(i,i)=strs(i,i)-extstrs(i)
     enddo   
     smax=0.d0
     do i=1,3
     do j=1,3
        if (i==j) then
           smax=smax+strs(i,j)**2
        else   
           smax=smax+strs(i,j)**2*2.d0
        endif   
     enddo   
     enddo   
     smax=dsqrt(smax)
! debug: >
     write(601,*)'strs'
     write(601,'("loopc, smax",i5,2e16.6,e23.10)')loopc,smax,sth,tote
     do i=1,3
        write(601,"(3f12.6)")(strs(i,j),j=1,3)
     enddo
! debug: <
     if (smax<sth) then
        write(*,'("strs converged. loopc, smax",i5,e12.4)')loopc,smax
        exit
     else
        vuv=vuv+0.5d0*tstep*matmul(strs,uv)/mcell        
     endif   
  enddo ! loopc

  call output_struc(901)
  close(901)
  write(6,*)'!!!! normaly end !!!!'

contains

subroutine show_argument()
  implicit none
  write(6,*)
  write(6,*)'argument:'
  write(6,*)'thisprogram inputfilename'
  write(6,*)
end subroutine show_argument
subroutine input
  use keyvalue

  integer :: ifin,ret,ilength,icolumn,iatompos
  integer :: i,j,itmp
  real(8) :: uvin(3,3),vtmp(3),aw
  character :: ctmp*2

  character(120):: inputfilename
  integer:: narg
 narg=command_argument_count()
  if (narg/=1) then 
     call show_argument()
     stop 100
  endif
  call get_command_argument(1,inputfilename)
  write(6,*)'input=',trim(inputfilename)
! unit cell
  call getkeyvalue(inputfilename,"unit_vec",unit=ifin,status=ret)
  write(*,*)'unit_vec'
  do i=1,3
     read(ifin,*)ctmp,(uvin(i,j),j=1,3) ! hm in QMAS
     write(*,"(3f23.16)")(uvin(i,j),j=1,3)
  enddo   
  close(ifin)

  call getkeyvalue(inputfilename,"unit_vec_column",icolumn,default=1)
  if (icolumn==2) then
     uv=transpose(uvin)
     write(*,*)'unit vector (column/raw): raw'
  else
     uv=uvin
     write(*,*)'unit vector (column/raw): column'
  endif   

  call getkeyvalue(inputfilename,"unit_length",ilength,default=1)
  if (ilength==1) uv=uv/bohr

! u1=uv(1:3,1)
! u2=uv(1:3,2)
! u3=uv(1:3,3)
  write(*,*)'unit vectors [Bohr]'
  write(*,"('u1 =',3f23.16)")uv(1:3,1)
  write(*,"('u2 =',3f23.16)")uv(1:3,2)
  write(*,"('u3 =',3f23.16)")uv(1:3,3)

  call cross_x(uv(:,1),uv(:,2),bv(:,3))
  call cross_x(uv(:,2),uv(:,3),bv(:,1))
  call cross_x(uv(:,3),uv(:,1),bv(:,2))
  omega=dot_product(bv(:,3),uv(:,3))
  write(*,*)'omega =',omega
  omegai=1.d0/omega
  bv=bv*omegai

! informtion of atoms
  call getkeyvalue(inputfilename,"number_atom",natom,default=0)
  write(*,*)'number of atoms =',natom
  allocate(ra(3,natom),rr(3,natom),rro(3,natom,2),iposfix(natom))
  allocate(frc(3,natom),vrr(3,natom))
  frc=0.d0
  vrr=0.d0

  call getkeyvalue(inputfilename,"number_element",nkatm,default=0)
  write(*,*)'number of elements =',nkatm
  allocate(zatm(natom),katm(natom),mass(natom),mfac(natom))
  mconv=(1.6605655d-27/9.109534d-31)*0.5d0

  call getkeyvalue(inputfilename,"atom_pos",iatompos,default=1)

  call getkeyvalue(inputfilename,"atom_list",unit=ifin,status=ret)
  do i=1,natom
     read(ifin,*)zatm(i),vtmp(1:3),katm(i),iposfix(i)
     if (iatompos==1) then
        rr(:,i)=vtmp
        ra(:,i)=matmul(uv,vtmp)
     else
        if (ilength==1) vtmp=vtmp/bohr
        ra(:,i)=vtmp
        rr(:,i)=matmul(vtmp,bv)
     endif   
  enddo   
  close(ifin)

  write(*,*)'atom_list: lattice coordinate'
  do i=1,natom
     write(*,"(i5,3f23.16,i4)")zatm(i),(rr(j,i),j=1,3),iposfix(i)
  enddo

  write(*,*)'atom_list: Cartesian coordinate [Bohr]'
  do i=1,natom
     write(*,"(i5,3f23.16,i4)")zatm(i),(ra(j,i),j=1,3),iposfix(i)
  enddo

! optimization
  call getkeyvalue(inputfilename,"number_max_relax",nloopa,default=10)
  write(*,*)'number_max_relax =',nloopa

  call getkeyvalue(inputfilename,"md_mode",imd,default=0)
  write(*,*)'md_mode =',imd

  call getkeyvalue(inputfilename,"number_max_relax_cell",nloopc,default=0)
  write(*,*)'number_max_relax_cell =',nloopc

  call getkeyvalue(inputfilename,"external_stress_v",extstrs,3,default=(/0.0d0,0.0d0,0.0d0/))
  write(*,"(a17,3f12.6)")'external_stress_v =',extstrs

  call getkeyvalue(inputfilename,"time_step",tstep,default=300.d0)
  write(*,*)'time_step =',tstep

  call getkeyvalue(inputfilename,"th_force",fth,default=5.d-5)
  write(*,*)'th_force =',fth

  call getkeyvalue(inputfilename,"th_stress",sth,default=5.d-7)
  write(*,*)'th_stress =',sth

! lattice_mass in QMAS
  call getkeyvalue(inputfilename,"mass_cell",mcell,default=5.d-4)
  write(*,*)'mass_cell =',mcell ! mltc in QMAS

! cell_opt_mode in QMAS
  call getkeyvalue(inputfilename,"md_mode_cell",imdc,default=2)
  write(*,*)'md_mode_cell =',imdc ! irlattice in QMAS

! atom mass
  if (imd<0) then
     do i=1,natom
        mass(i)=12.d0*mconv
        mfac(i)=1.d0
     enddo   
  else   
     do i=1,natom
        call massset(zatm(i),aw)
        mass(i)=aw*mconv
        mfac(i)=1.d0
     enddo   
  endif

! force field
  call getkeyvalue(inputfilename,"force_field",ifrcf,default=0)
  write(*,*)'force_field =',ifrcf

end subroutine input      

subroutine updt_cell

  integer :: ifin,ret,ilength,icolumn,iatompos

  uvo(:,:,2)=uvo(:,:,1)
  uvo(:,:,1)=uv
  strso(:,:,2)=strso(:,:,1)
  strso(:,:,1)=strs

  vuv=vuv+0.5d0*tstep*matmul(strs,uv)/mcell

  if (imdc==1) then ! steepest descent
     call latticerelax_sd
  elseif (imdc==2) then ! quenched MD
     call latticerelax
  else   
     call lattice_simple_relax
  endif   

  call cross_x(uv(:,1),uv(:,2),bv(:,3))
  call cross_x(uv(:,2),uv(:,3),bv(:,1))
  call cross_x(uv(:,3),uv(:,1),bv(:,2))
  omega=dot_product(bv(:,3),uv(:,3))
  write(*,*)'omega =',omega
  omegai=1.d0/omega
  bv=bv*omegai

end subroutine updt_cell    

subroutine tote_frc_strs
  real*8 faemp

  if (ifrcf==1) then ! Stillinger-Weber for Si
     call Stillinger_Weber
  else
     stop 'etot_frc_strs error'
  endif

  fmax=0.d0
  do i=1,natom
     faemp=dsqrt(sum(frc(:,i)**2))
     if (iposfix(i)/=0.and.faemp>fmax) then
        fmax=faemp
        iamax=i
     endif
  enddo

  write(*,*)'tote,fmax =',tote,fmax

end subroutine tote_frc_strs

subroutine Stillinger_Weber

  integer :: i,j,k
  real*8 :: tote2,tote3,frc2(3,natom),frc3(3,natom),strs2(3,3),strs3(3,3)
  real*8 :: rri(3),rrj(3),rrk(3)
  real*8 :: rrij(3),raij(3),dij,rrik(3),raik(3),dik,rrjk(3),rajk(3),djk
  real*8 :: toteij,frcijr,toteijk,frcijk(3,3),vtmp(3)
! debug:
  real*8 :: sigma,eps,e1,e2,e3,vol,rho,faemp
  parameter (sigma=3.959164919d0) ! in Bohr
  parameter (eps=7.968005097d-2) ! In Hartree

! two-body term
  tote2=0.d0
  frc2=0.d0
  strs2=0.d0
!$omp parallel do private(rri,rrj,rrij,raij,dij,frcijr,vtmp,frc2), reduction(+:tote2,strs2) ,  schedule(dynamic,2)
  do i=1,natom-1
  do j=i+1,natom
     rri=rr(:,i)
     rrj=rr(:,j)
     call frac_diff_min(rri,rrj,rrij)
     raij=matmul(uv,rrij)
     dij=sqrt(sum(raij**2))
     call SW_Si_2body(dij,toteij,frcijr)
     tote2=tote2+toteij
     vtmp= frcijr*raij(:)/dij
     frc2(:,i)= frc2(:,i)+vtmp
     frc2(:,j)= frc2(:,j)-vtmp
     do k=1,3
     strs2(:,k)=strs2(:,k)+vtmp*raij(k)
     enddo
  enddo ! j
  enddo ! i

! three-body term
  tote3=0.d0
  frc3=0.d0
  strs3=0.d0
!$omp parallel do private(rri,rrj,rrij,raij,rrik,raik,rrk,rrjk,rajk,toteijk,frcijk) , reduction(+:tote3,frc3) , schedule(dynamic,2)
  do i=1,natom-2
  do j=i+1,natom-1
  do k=j+1,natom
     rri=rr(:,i)
     rrj=rr(:,j)
     call frac_diff_min(rri,rrj,rrij)
     raij=matmul(uv,rrij)
     rri=rr(:,i)
     rrk=rr(:,k)
     call frac_diff_min(rri,rrk,rrik)
     raik=matmul(uv,rrik)
     rrj=rr(:,j)
     rrk=rr(:,k)
     call frac_diff_min(rrj,rrk,rrjk)
     rajk=matmul(uv,rrjk)
     call SW_Si_3body(raij,raik,rajk,toteijk,frcijk,strs3)
     tote3=tote3+toteijk
     frc3(:,i)= frc3(:,i)+frcijk(:,1)
     frc3(:,j)= frc3(:,j)+frcijk(:,2)
     frc3(:,k)= frc3(:,k)+frcijk(:,3)
  enddo ! k
  enddo ! j
  enddo ! i

! total
  tote=tote2+tote3
  frc=frc2+frc3
! debug
!  strs=strs2+strs3
  strs=strs2-strs3

! debug:
!  e1=tote/dble(natom)/eps
!  e2=tote2/dble(natom)/eps
!  e3=tote3/dble(natom)/eps
!  vol=omega/dble(natom)
!  rho=1.d0/vol*sigma**3
!  write(902,"(2i5,2x,4f12.6)")loopc,loopa,rho,e1,e2,e3
!  write(903,*)'frc'
!  write(903,*)'*** loopa,loopc',loopa,loopc
!  do i=1,natom
!     write(903,"(i5,3f12.6)")i,frc(1:3,i)
!  enddo   
!
!  fmax=0.d0
!  do i=1,natom
!     faemp=dsqrt(sum(frc(:,i)**2))
!     if (iposfix(i)/=0.and.faemp>fmax) then
!        fmax=faemp
!        iamax=i
!     endif
!  enddo
!  write(906,*)'loopa, tote, fmax =',loopa,e1,fmax
! debug: <

end subroutine Stillinger_Weber

! two body contributions to Si-Si interaction 
! Stillinger and Weber, Phys.Rev.B31,5262(1985)
! atomic unit
subroutine SW_Si_2body(r,tote2,frc2r)
  implicit none
  real*8 A,B
  parameter (A=7.049556277d0)
  parameter (B=0.6022245584d0)
  integer p,q
  parameter (p=4)
  parameter (q=0)
  real*8 a_cut,sigma,eps
  parameter (a_cut=1.80d0)
  parameter (sigma=3.959164919d0) ! in Bohr
  parameter (eps=7.968005097d-2) ! in Hartree
  real*8 r,tote2,frc2r
  real*8 r0,f2,df2  
!=======================================================================
  r0=r/sigma
  if (r0.ge.a_cut) then 
     tote2=0d0
     frc2r=0d0
  else 
     f2 = A*(B*r0**(-p)-r0**(-q))*dexp(1d0/(r0-a_cut))
     df2 = A*(-p*B*r0**(-p-1)+q*r0**(-q-1))*dexp(1d0/(r0-a_cut)) &
          -A*(B*r0**(-p)-r0**(-q))*dexp(1d0/(r0-a_cut))*(r0-a_cut)**(-2)
     tote2 = eps*f2
     frc2r = -eps*df2/sigma
  endif
  return
end subroutine SW_Si_2body

! three body contributions to Si-Si-Si interaction 
! Stillinger and Weber, Phys.Rev.B31,5262(1985)
! atomic unit
subroutine SW_Si_3body(r12,r13,r23,tote3,frc3,strs3)
  implicit none
  real*8 a_cut,sigma,eps
  parameter (a_cut=1.80d0)
  parameter (sigma=3.959164919d0) ! in Bohr
  parameter (eps=7.968005097d-2) ! in Hartree
  real*8 r12(3),r13(3),r23(3),s12(3),s13(3),s23(3),v1(3),v2(3)
  real*8 tote3,frc3(3,3),strs3(3,3)
  real*8 d12,d13,d23
  real*8 h,dh(3,3)
  integer i,j,k
!=======================================================================
  s12=r12/sigma
  s13=r13/sigma
  s23=r23/sigma

  d12=dsqrt(sum(s12**2))
  d13=dsqrt(sum(s13**2))
  d23=dsqrt(sum(s23**2))

  tote3=0.d0
  frc3=0.d0
   
  if (d12<a_cut.and.d13<a_cut) then
     v1=-s12
     v2=-s13
     call sw3body(a_cut,d12,d13,v1,v2,h,dh)
!     call sw3body(a_cut,d12,d13,-s12,-s13,h,dh)
     tote3=tote3+eps*h
     frc3(:,1)=frc3(:,1)-eps/sigma*dh(:,1)
     frc3(:,2)=frc3(:,2)-eps/sigma*dh(:,2)
     frc3(:,3)=frc3(:,3)-eps/sigma*dh(:,3)
     do k=1,3
     strs3(:,k)=strs3(:,k)-eps/sigma*dh(:,2)*r12(k)-eps/sigma*dh(:,3)*r13(k)
     enddo
  endif   

  if (d12<a_cut.and.d23<a_cut) then
     v1=s12
     v2=-s23
     call sw3body(a_cut,d12,d23,v1,v2,h,dh)
!     call sw3body(a_cut,d12,d23,s12,-s23,h,dh)
     tote3=tote3+eps*h
     frc3(:,2)=frc3(:,2)-eps/sigma*dh(:,1)
     frc3(:,1)=frc3(:,1)-eps/sigma*dh(:,2)
     frc3(:,3)=frc3(:,3)-eps/sigma*dh(:,3)
     do k=1,3
     strs3(:,k)=strs3(:,k)+eps/sigma*dh(:,2)*r12(k)-eps/sigma*dh(:,3)*r23(k)
     enddo
  endif   

  if (d13<a_cut.and.d23<a_cut) then
     v1=s13
     v2=s23
     call sw3body(a_cut,d13,d23,v1,v2,h,dh)
!     call sw3body(a_cut,d13,d23,s13,s23,h,dh)
     tote3=tote3+eps*h
     frc3(:,3)=frc3(:,3)-eps/sigma*dh(:,1)
     frc3(:,1)=frc3(:,1)-eps/sigma*dh(:,2)
     frc3(:,2)=frc3(:,2)-eps/sigma*dh(:,3)
     do k=1,3
     strs3(:,k)=strs3(:,k)+eps/sigma*dh(:,2)*r13(k)+eps/sigma*dh(:,3)*r23(k)
     enddo
  endif   

  return
  end subroutine SW_Si_3body
!***********************************************************************
! h and its derivatives in Stillinger-Weber potential
subroutine sw3body(a_cut,d12,d13,s12,s13,h,dh)
  implicit none

  real*8 lamda,gamma
  parameter (lamda=21.0d0)
  parameter (gamma=1.20d0)
  real*8 a_cut
  real*8 d12,d13,s12(3),s13(3)
  real*8 h,dh(3,3)
  real*8 b,c
!=======================================================================
  b = sum(s12*s13)
  c = b/(d12*d13) + 1d0/3d0

  if (c==0d0) then 
     h = 0d0
     dh = 0.d0
  else
     h = lamda*dexp(gamma/(d12-a_cut)+gamma/(d13-a_cut))*c**2
     dh(:,1) = h*(gamma*s12/d12/(d12-a_cut)**2 & 
                + gamma*s13/d13/(d13-a_cut)**2) &        
             + 2d0*h/c*(-(s12+s13)/(d12*d13)-b*(d12*(-s13/d13) &
             + d13*(-s12/d12))/(d12*d13)**2)
     dh(:,2) = h*(gamma*(-s12/d12)/(d12-a_cut)**2) &
             + 2d0*h/c*(s13/(d12*d13)-b*(s12/d12)/(d12**2*d13))
     dh(:,3) = h*(gamma*(-s13/d13)/(d13-a_cut)**2) &
             + 2d0*h/c*(s12/(d12*d13)-b*(s13/d13)/(d12*d13**2))
  endif

  return
end subroutine sw3body
!***********************************************************************
#if 1
#include "frac_diff_min.f90"
#else
subroutine frac_diff_min(p,p1,diff)
   implicit none
   real(8),intent(in):: p(3)
   real(8),intent(out)::p1(3)
   integer:: i1,i2,i3,ilist(3)
   real(8):: diffmin(3),diff(3),difflen

   diff=p-p1
   difflen =  sum(diff**2)
   ilist=[0,0,0]
   do i3=-1,1
   do i2=-1,1
   do i1=-1,1
     diff = p-p1-[ i1,i2,i3] 
     if (difflen > sum(diff**2) ) then 
        difflen= sum(diff**2) 
        ilist= [ i1,i2,i3]
     endif
   enddo
   enddo
   enddo
     
   i1=ilist(1); i2=ilist(2); i3=ilist(3)
   p1= p1+ [ real(i1,kind=8),real(i2,kind=8),real(i3,kind=8) ]
   diff= p-p1

end subroutine frac_diff_min
#endif

subroutine latticerelax

  real*8 :: auv(3,3,2),dtmp,c1,c2
  integer :: i,j

  do i=1,2
    auv(:,:,i)=matmul(strso(:,:,i),uvo(:,:,i))
  enddo

  do j=1,3
  do i=1,3
     if ((auv(i,j,1)*auv(i,j,2))<(-tol*1.d-10)) then
        if (loopc>1) then
!           dtmp=-auv(i,j,2)*uvo(i,j,1)/(auv(i,j,1)-auv(i,j,2)) &
!                +auv(i,j,1)*uvo(i,j,2)/(auv(i,j,1)-auv(i,j,2)) &
!                -uv(i,j)
           c1=-auv(i,j,2)/(auv(i,j,1)-auv(i,j,2))
           c1=dmin1(c1,1.d0)
           c1=dmax1(c1,0.d0)
           c2=1.d0-c1
           dtmp=c1*uvo(i,j,1)+c2*uvo(i,j,2)-uv(i,j)
        else   
           dtmp=0.5d0*tstep*tstep*auv(i,j,1)/mcell
        endif
        dtmp=dmax1(dtmp,-0.1d0)
        dtmp=dmin1(dtmp, 0.1d0)
        uv(i,j)=uv(i,j)+dtmp
        vuv(i,j)=0.d0
     else
        uv(i,j)=uv(i,j)+tstep*vuv(i,j)
     endif   
  enddo ! i   
  enddo ! j   

end subroutine latticerelax

subroutine lattice_simple_relax

 uv=uv+tstep*vuv
 vuv=0.d0

end subroutine lattice_simple_relax

subroutine latticerelax_sd
! modified from latticerelax_cg in QMAS
  integer i,j,k
  real*8 lambt1,lambt2,lambu1,lambu2,lambda,dtmp1,dtmp2,gamma
  real*8 acoe,bcoe,dhmmax,lamb0,dd

!  trial step
  lambt1=25.0d0
  lambt2=50.0d0

!  maximum of lambda
  lambu1=100.0d0
  lambu2=10.0d0

  if(loopc==2) then
     isd_cell=0
     guv(:,:,0)=0.0d0
  endif

  guv(:,:,2)=matmul(strs,uv)

  if(loopc>=3) then
     dtmp1=0.0d0
     dtmp2=0.0d0
     do j=1,3
     do i=1,3
        dtmp1=dtmp1+guv(i,j,0)*guv(i,j,0)
        dtmp2=dtmp2+guv(i,j,0)*guv(i,j,2)
     end do
     end do
     write(6,'(1x,"g0*s_cell=",d16.8,2x,"g2*s_cell=",d16.8)') dtmp1,dtmp2
     if(dabs(dtmp2).lt.dabs(dtmp1)*0.01d0) isd_cell=0
  endif

  if(isd_cell==4) isd_cell=0

  isd_cell=isd_cell+1

! --- new search direction (sdv)

  if(isd_cell==1) then

     if(loopc<=3) then
        lambda=lambt1
     else
        lambda=lambt2
     endif

     dd=0.0d0
     do j=1,3
     do i=1,3
        dd=dmax1(dd,dabs(guv(i,j,2)))
     enddo
     enddo
     lambda=dmin1(lambda,(0.5d0/(dd*bohr)))

     uv(:,:)=uv(:,:)+lambda*guv(:,:,2)
     duv0(:,:)=lambda*guv(:,:,2)

     guv(:,:,0)=guv(:,:,2)
     guv(:,:,1)=guv(:,:,2)

! --- search the minimum along the search direction (sdv)

  else if(isd_cell >= 2) then

     dtmp1=0.0d0
     dtmp2=0.0d0
     do j=1,3
     do i=1,3
        dtmp1=dtmp1+guv(i,j,0)*guv(i,j,1)
        dtmp2=dtmp2+guv(i,j,0)*guv(i,j,2)
     enddo
     enddo

     if(dtmp1 < 0.0d0) then
        dtmp1=-dtmp1
        dtmp2=-dtmp2
     endif

     write(6,'(1x,"g1*s_cell=",d16.8,2x,"g2*s_cell=",d16.8)')dtmp1,dtmp2

     if(dtmp1 > dtmp2) then
        acoe=-dtmp1
        bcoe=0.5d0*(-dtmp2-acoe)
        lambda=-0.5d0*acoe/bcoe
        if(lambda > lambu2) then
           lambda=lambu2-1.0d0
        elseif(lambda < (-lambu2)) then
           lambda=-lambu2-1.0d0
        else
           lambda=lambda-1.0d0
        endif
     else
        lambda=50.0d0-1.0d0
     endif

     uv(:,:)=uv(:,:)+lambda*duv0(:,:)
     duv0(:,:)=lambda*duv0(:,:)

     guv(:,:,1)=guv(:,:,2)

  end if ! isd_cell

  write(6,'(1x,"isd_cell=",i2,2x,"lambda=",d16.8)') isd_cell,lambda

  write(6,*)

  return
end subroutine latticerelax_sd

subroutine updt_coord

  integer :: ia
  real*8 :: rfac,p1,p2,vrrmax,vrrabs
  real*8,allocatable :: ratmp(:,:,:)

  if (loopa>1)rro(:,:,2)=rro(:,:,1)
  rro(:,:,1)=rr(:,:)
  vrrmax=0.d0
  if (loopa>1) then
     do ia=1,natom
!        rfac=0.5d0*tstep*omegai/(mass(ia)*mfac(ia))
        rfac=0.5d0*tstep/(mass(ia)*mfac(ia))
        vrr(:,ia)=vrr(:,ia)+rfac*matmul(frc(:,ia),bv(:,:))
        vrrabs=dsqrt(sum(vrr(:,ia)**2))
        if (vrrabs>vrrmax) vrrmax=vrrabs
     enddo   
  endif   

! debug:
!  write(*,*)'vrrmax =',loopa,vrrmax

  call atomrelax

  if (loopa>2) then
     allocate(ratmp(3,natom,0:2))
     do ia=1,natom
        ratmp(:,:,0)=matmul(uv,rr)
        ratmp(:,:,1)=matmul(uv,rro(:,:,1))
        ratmp(:,:,2)=matmul(uv,rro(:,:,2))
     enddo   
     p1=0.d0
     p2=0.d0
     do ia=1,natom
        p1=p1+sum((ratmp(:,ia,0)-ratmp(:,ia,1))*(ratmp(:,ia,1)-ratmp(:,ia,2)))
        p2=p2+sum((ratmp(:,ia,1)-ratmp(:,ia,2))**2)
     enddo   
     deallocate(ratmp)
     if (dabs(p2)<1.d-8) then
        alphalm=0.d0
     else
        alphalm=p1/p2
     endif
  else ! loopa
     alphalm=0.d0
  endif ! loopa>2   

end subroutine updt_coord    

subroutine output_struc(ifo)
  integer :: ifo,i,j

  write(ifo,*)'*** unit vectors [Bohr]: loopc =',loopc
  do i=1,3
     write(ifo,"(3f23.16)")(uv(i,j),j=1,3)
  enddo
  write(ifo,*)'*** internal lattice coordinates: loopc =',loopc
  do j=1,natom
     write(ifo,"(3f23.16)")(rr(i,j),i=1,3)
  enddo

  return
end subroutine output_struc

end program main
!--------1---------2---------3---------4---------5---------6---------7--
      subroutine a2b(a,b)
      implicit none
      real*8 :: pi,a(3,3),b(3,3),ab
! sum[i] a(i,j)*b(i,k) = delta(j,k)

!      pi = 4d0*datan(1.0d0)
      call cross_x(a(1,2),a(1,3),b(1,1))
      call cross_x(a(1,3),a(1,1),b(1,2))
      call cross_x(a(1,1),a(1,2),b(1,3))
      ab = sum(a(1:3,1)*b(1:3,1))

!      b(:,1) = 2d0*pi*b(:,1)/ab
      b(:,1) = b(:,1)/ab
      ab = sum(a(1:3,2)*b(1:3,2))
!      b(:,2) = 2d0*pi*b(:,2)/ab
      b(:,2) = b(:,2)/ab
      ab = sum(a(1:3,3)*b(1:3,3))
!      b(:,3) = 2d0*pi*b(:,3)/ab
      b(:,3) = b(:,3)/ab

      return
      end
!--------1---------2---------3---------4---------5---------6---------7--
      subroutine cross_x(a,b,c)
      implicit none
      real(8) :: a(3),b(3),c(3)
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)
      return
      end
!--------1---------2---------3---------4---------5---------6---------7--
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
