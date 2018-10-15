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

module frcfield
  use paramlist
  implicit none

  contains

!--------1---------2---------3---------4---------5---------6---------7--
!--------1---------2---------3---------4---------5---------6---------7--
subroutine Stillinger_Weber
  use periodic_lattice_module, only: cell_list_type, kr, periodic_lattice_type

  integer :: i,j,k
  integer :: jj,kk
  real*8 :: tote2,tote3,frc2(3,QMD%natom),frc3(3,QMD%natom),strs2(3,3),strs3(3,3)
  real*8 :: rri(3),rrj(3),rrk(3)
  real*8 :: rrij(3),raij(3),dij,rrik(3),raik(3),rrjk(3),rajk(3)
  real*8 :: toteij,frcijr,toteijk,frcijk(3,3),vtmp(3),strsijk(3,3)
  type(periodic_lattice_type) :: lattice
  type(cell_list_type), allocatable :: cell_of_replica(:,:)
! debug:
  real*8 :: a_cut
  real*8 :: sigma,eps,e1,e2,e3,vol,rho,faemp
  parameter (a_cut=1.8d0)
  parameter (sigma=3.959164919d0) ! in Bohr
  parameter (eps=7.968005097d-2) ! in Hartree

  lattice%direct_lattice=real(QMD%uv,kr)
  allocate(cell_of_replica(QMD%natom,QMD%natom))

! two-body term
!$omp parallel do private(j)
  do i=1,QMD%natom
     do j=1,QMD%natom
        cell_of_replica(j,i)=lattice%get_cell_of_replica(real(QMD%rr(:,j),kr), real(QMD%rr(:,i),kr), real(a_cut*sigma,kr))
     enddo
  enddo
  tote2=0.d0
  frc2=0.d0
  strs2=0.d0
!$omp parallel do private(rri,rrj,rrij,raij,dij,toteij,frcijr,vtmp,j,jj,k), reduction(+:tote2,frc2,strs2) ,  schedule(dynamic,2)
  do i=1,QMD%natom
     rri=QMD%rr(:,i)
     do j=i,QMD%natom
        do jj=1,size(cell_of_replica(j,i)%list,dim=2)
           if (j==i.and.all(cell_of_replica(j,i)%list(:,jj)==floor(rri(:)))) cycle
           rrj=cell_of_replica(j,i)%list(:,jj)+modulo(QMD%rr(:,j),1.0d0)
           rrij=rrj(:)-rri(:)
           raij=matmul(QMD%uv,rrij)
           dij=sqrt(sum(raij**2))
           if (QMD%zatm(i)==14.and.QMD%zatm(j)==14) &
           call SW_Si_2body(dij,toteij,frcijr)
           if (i==j) then
              toteij=toteij/2.0d0
              frcijr=frcijr/2.0d0
           endif
           tote2=tote2+toteij
           vtmp= frcijr*raij(:)/dij
           frc2(:,i)= frc2(:,i)-vtmp
           frc2(:,j)= frc2(:,j)+vtmp
           do k=1,3
              strs2(:,k)=strs2(:,k)+vtmp*raij(k)
           enddo
        enddo ! jj
     enddo ! j
  enddo ! i

! three-body term
!$omp parallel do private(j)
  do i=1,QMD%natom
     do j=1,QMD%natom
        cell_of_replica(j,i)=lattice%get_cell_of_replica(real(QMD%rr(:,j),kr), real(QMD%rr(:,i),kr), real(2*a_cut*sigma,kr)) ! interaction range is twice the cutoff radius
     enddo
  enddo
  tote3=0.d0
  frc3=0.d0
  strs3=0.d0
!$omp parallel do private(rri,rrj,rrij,raij,rrik,raik,rrk,rrjk,rajk,toteijk,frcijk,strsijk,j,jj,k,kk) , reduction(+:tote3,frc3,strs3) , schedule(dynamic,2)
  do i=1,QMD%natom
     rri=QMD%rr(:,i)
     do j=i,QMD%natom
        do jj=1,size(cell_of_replica(j,i)%list,dim=2)
           if (j==i.and.all(cell_of_replica(j,i)%list(:,jj)==floor(rri(:)))) cycle
           rrj=cell_of_replica(j,i)%list(:,jj)+modulo(QMD%rr(:,j),1.0d0)
           rrij=rrj(:)-rri(:)
           raij=matmul(QMD%uv,rrij)
           do k=j,QMD%natom
              do kk=1,size(cell_of_replica(k,i)%list,dim=2)
                 if (k==i.and.all(cell_of_replica(k,i)%list(:,kk)==floor(rri(:)))) cycle
                 if (k==j.and.kk==jj) cycle
                 rrk=cell_of_replica(k,i)%list(:,kk)+modulo(QMD%rr(:,k),1.0d0)
                 rrik=rrk(:)-rri(:)
                 raik=matmul(QMD%uv,rrik)
                 rrjk=rrk(:)-rrj(:)
                 rajk=matmul(QMD%uv,rrjk)
                 if (QMD%zatm(i)==14.and.QMD%zatm(j)==14.and.QMD%zatm(k)==14) &
                 call SW_Si_3body(raij,raik,rajk,toteijk,frcijk,strsijk)
                 if ((i==j.and.k/=i).or.(j==k.and.i/=j).or.(k==i.and.j/=k)) then
                    toteijk=toteijk/2.0d0
                    frcijk=frcijk/2.0d0
                    strsijk=strsijk/2.0d0
                 elseif (i==j.and.j==k) then
                    toteijk=toteijk/6.0d0
                    frcijk=frcijk/6.0d0
                    strsijk=strsijk/6.0d0
                 endif
                 tote3=tote3+toteijk
                 frc3(:,i)= frc3(:,i)+frcijk(:,1)
                 frc3(:,j)= frc3(:,j)+frcijk(:,2)
                 frc3(:,k)= frc3(:,k)+frcijk(:,3)
                 strs3=strs3+strsijk
              enddo ! kk
           enddo ! k
        enddo ! jj
     enddo ! j
  enddo ! i

  deallocate(cell_of_replica)

! total
  QMD%tote=tote2+tote3
  QMD%frc=frc2+frc3
  QMD%strs=(strs2+strs3)/QMD%omega
!  QMD%strs=(strs2-strs3)/QMD%omega

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
  strs3=0.d0
   
  if (d12<a_cut.and.d13<a_cut) then
!     v1=-s12
!     v2=-s13
     v1=s12
     v2=s13
     call sw3body(a_cut,d12,d13,v1,v2,h,dh)
!     call sw3body(a_cut,d12,d13,-s12,-s13,h,dh)
     tote3=tote3+eps*h
     frc3(:,1)=frc3(:,1)-eps/sigma*dh(:,1)
     frc3(:,2)=frc3(:,2)-eps/sigma*dh(:,2)
     frc3(:,3)=frc3(:,3)-eps/sigma*dh(:,3)
     do k=1,3
      strs3(:,k)=strs3(:,k)-eps/sigma*dh(:,2)*r12(k)-eps/sigma*dh(:,3)*r13(k)
!     strs3(:,k)=strs3(:,k)+eps/sigma*dh(:,2)*r12(k)+eps/sigma*dh(:,3)*r13(k)
     enddo
  endif   

  if (d12<a_cut.and.d23<a_cut) then
!     v1=s12
!     v2=-s23
     v1=-s12
     v2=s23
     call sw3body(a_cut,d12,d23,v1,v2,h,dh)
!     call sw3body(a_cut,d12,d23,s12,-s23,h,dh)
     tote3=tote3+eps*h
     frc3(:,2)=frc3(:,2)-eps/sigma*dh(:,1)
     frc3(:,1)=frc3(:,1)-eps/sigma*dh(:,2)
     frc3(:,3)=frc3(:,3)-eps/sigma*dh(:,3)
     do k=1,3
      strs3(:,k)=strs3(:,k)+eps/sigma*dh(:,2)*r12(k)-eps/sigma*dh(:,3)*r23(k)
!     strs3(:,k)=strs3(:,k)-eps/sigma*dh(:,2)*r12(k)+eps/sigma*dh(:,3)*r23(k)
     enddo
  endif   

  if (d13<a_cut.and.d23<a_cut) then
!     v1=s13
!     v2=s23
     v1=-s13
     v2=-s23
     call sw3body(a_cut,d13,d23,v1,v2,h,dh)
!     call sw3body(a_cut,d13,d23,s13,s23,h,dh)
     tote3=tote3+eps*h
     frc3(:,3)=frc3(:,3)-eps/sigma*dh(:,1)
     frc3(:,1)=frc3(:,1)-eps/sigma*dh(:,2)
     frc3(:,2)=frc3(:,2)-eps/sigma*dh(:,3)
     do k=1,3
      strs3(:,k)=strs3(:,k)+eps/sigma*dh(:,2)*r13(k)+eps/sigma*dh(:,3)*r23(k)
!     strs3(:,k)=strs3(:,k)-eps/sigma*dh(:,2)*r13(k)-eps/sigma*dh(:,3)*r23(k)
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
!--------1---------2---------3---------4---------5---------6---------7--
!--------1---------2---------3---------4---------5---------6---------7--
subroutine Tsuneyuki

  integer :: i,j,k
  real*8 :: tote2,frc2(3,QMD%natom),strs2(3,3)
  real*8 :: toteEwald,frcEwald(3,QMD%natom),strsEwald(3,3)
  real*8 :: rri(3),rrj(3)
  real*8 :: rrij(3),raij(3),dij
  real*8 :: toteij,frcijr,vtmp(3)
  real*8 :: toteijEwald,frcijEwald(3),strsijEwald(3,3)
  !real*8 :: sigma,eps,e1,e2,e3,vol,rho,faemp
  !parameter (sigma=3.959164919d0) ! in Bohr
  !parameter (eps=7.968005097d-2) ! In Hartree

  tote2=0.d0
  frc2=0.d0
  strs2=0.d0
!$omp parallel do private(rri,rrj,rrij,raij,dij,frcijr,vtmp,j,k), reduction(+:tote2,frc2,strs2) ,  schedule(dynamic,2)
  do i=1,QMD%natom-1
     if (QMD%zatm(i)/=8.and.QMD%zatm(i)/=14) stop 'Tsuneyuki: zatm(i) error'
  do j=i+1,QMD%natom
     if (QMD%zatm(j)/=8.and.QMD%zatm(j)/=14) stop 'Tsuneyuki: zatm(j) error'
     if (QMD%zatm(i)==14.and.QMD%zatm(j)==14) cycle
     rri=QMD%rr(:,i)
     rrj=QMD%rr(:,j)
!     call frac_diff_min(rri,rrj,rrij) ! rrij = rri - rrj
     call frac_diff_min(rrj,rri,rrij) ! rrij = rrj - rri
     raij=matmul(QMD%uv,rrij)
     dij=sqrt(sum(raij**2))
     call Tsuneyuki_Uij(QMD%zatm(i),QMD%zatm(j),dij,toteij,frcijr)
     tote2=tote2+toteij
     vtmp= frcijr*raij(:)/dij
     frc2(:,i)= frc2(:,i)-vtmp
     frc2(:,j)= frc2(:,j)+vtmp
     do k=1,3
     strs2(:,k)=strs2(:,k)+vtmp*raij(k)
     enddo
  enddo ! j
  enddo ! i

  toteEwald=0.d0
  frcEwald=0.d0
  strsEwald=0.d0
!$omp parallel do private(rri,rrj,rrij,raij,toteijEwald,frcijEwald,strsijEwald,j), reduction(+:toteEwald,frcEwald,strsEwald) ,  schedule(dynamic,2)
  do i=1,QMD%natom
     if (QMD%zatm(i)/=8.and.QMD%zatm(i)/=14) stop 'Tsuneyuki: zatm(i) error'
  do j=1,QMD%natom
     if (QMD%zatm(j)/=8.and.QMD%zatm(j)/=14) stop 'Tsuneyuki: zatm(j) error'
     if (QMD%zatm(i)==14.and.QMD%zatm(j)==14) cycle
     rri=QMD%rr(:,i)
     rrj=QMD%rr(:,j)
     rrij=rrj-rri
     raij=matmul(QMD%uv,rrij)
     call Tsuneyuki_Ewald(QMD%zatm(i),QMD%zatm(j),raij,toteijEwald,frcijEwald,strsijEwald)
     toteEwald=toteEwald+0.5d0*toteijEwald
     frcEwald(:,i)=frcEwald(:,i)-frcijEwald
     strsEwald=strsEwald+0.5d0*strsijEwald
  enddo ! j
  enddo ! i

! total
  QMD%tote=tote2+toteEwald
  QMD%frc=frc2+frcEwald
  QMD%strs=(strs2+strsEwald)/QMD%omega

!  write(*,*)'frc'
!  do i=1,QMD%natom
!     write(*,'(i5,3f12.6)')i,QMD%frc(1:3,i)
!  enddo

end subroutine Tsuneyuki

! Tsuneyuki et al., Phys.Rev.Lett.61,869(1988)
! atomic unit
subroutine Tsuneyuki_Uij(zatmi,zatmj,r,tote2,frc2r)
  implicit none
  real*8 f0,etainv,eta
  parameter (f0=1.0d0)
  parameter (etainv=1.4d0)
  real*8 rminSiO,rminOO
  parameter (rminSiO=1.3d0, rminOO=1.0d0)
  integer zatmi,zatmj
  real*8 r,tote2,frc2r
  logical:: lpair
!
  integer zi,zj
  real*8 ai,aj,bi,bj,ci,cj,qi,qti,qj,qtj,UijCoul,dUijCoul,gij,dgij
  real*8 Uij,dUij,Uij2,dUij2
  real*8 rang,rfac
!=======================================================================
  rang=r*bohr
  lpair=.false.
!  if (r.ge.a_cut) then 
!     tote2=0d0
!     frc2r=0d0
!  else 
     call Tsuneyuki_Param(zatmi,ai,bi,ci,qi,qti)
     call Tsuneyuki_Param(zatmj,aj,bj,cj,qj,qtj)
     eta=bohr/etainv
     if ((zatmi.eq.8.and.zatmj.eq.14).or.(zatmi.eq.14.and.zatmj.eq.8)) then
        call Tsuneyuki_gSiO(r,eta,gij,dgij)
        lpair=.true.
     endif  
     if (zatmi.eq.8.and.zatmj.eq.8) then
        call Tsuneyuki_gOO(r,eta,gij,dgij)
        lpair=.true.
     endif  

!     if (lpair==.false.) then
     if (.not.lpair) then
        tote2=0.d0
        frc2r=0.d0
        return
     endif   

! Non Coulomb part: in units of Angstrom and kcal/mol
     rfac=1.d0
     if ((zatmi.eq.8.and.zatmj.eq.14).or.(zatmi.eq.14.and.zatmj.eq.8)) then
        if (rang<rminSiO) rfac=0.0d0
     elseif (zatmi.eq.8.and.zatmj.eq.8) then
        if (rang<rminOO) rfac=0.0d0
     endif   
     Uij2=f0*(bi+bj)*exp((ai+aj-rang)/(bi+bj))-ci*cj/(rang**6)*rfac
     dUij2=-f0*(bi+bj)*exp((ai+aj-rang)/(bi+bj))/(bi+bj)+6.d0*ci*cj/(rang**7)*rfac

! Non Coulomb part: convert to atomic units
     Uij2=Uij2/kcalpermol
     dUij2=dUij2*bohr/kcalpermol

! Coulomb(short) part: atomic units
!     UijCoul=qti*qtj*(1.d0-gij)/r+qi*qj*gij/r
!     dUijCoul=-qti*qtj*(1.d0-gij)/(r**2)-qti*qtj*dgij/r &
!              -qi*qj*gij/(r**2)+qi*qj*dgij/r
     UijCoul=-qti*qtj*gij/r+qi*qj*gij/r
     dUijCoul=qti*qtj*gij/(r**2)-qti*qtj*dgij/r-qi*qj*gij/(r**2)+qi*qj*dgij/r
!
     Uij=UijCoul+Uij2
!!     Uij=UijCoul
!!     Uij=Uij2
     dUij=dUijCoul+dUij2
     tote2=Uij
     frc2r=-dUij
!  endif
  return
end subroutine Tsuneyuki_Uij

subroutine Tsuneyuki_Ewald(zatmi,zatmj,rij,toteEwald,frcEwald,strsEwald)
  use ewald_module, only: M_PI, calcPotentialEwald, calcForceEwald, calcStressEwald
  implicit none
  integer,intent(in):: zatmi,zatmj
  real*8,intent(in):: rij(3)
  real*8,intent(out):: toteEwald,frcEwald(3),strsEwald(3,3)
  logical:: lpair
!
  real*8 ai,aj,bi,bj,ci,cj,qi,qti,qj,qtj
!=======================================================================
  lpair=.false.
  if ((zatmi.eq.8.and.zatmj.eq.14).or.(zatmi.eq.14.and.zatmj.eq.8)) then
    lpair=.true.
  endif
  if (zatmi.eq.8.and.zatmj.eq.8) then
    lpair=.true.
  endif

  if (.not.lpair) then
    toteEwald=0.d0
    frcEwald=0.d0
    return
  endif

  call Tsuneyuki_Param(zatmi,ai,bi,ci,qi,qti)
  call Tsuneyuki_Param(zatmj,aj,bj,cj,qj,qtj)

! Coulomb(long) part: atomic units
  call calcPotentialEwald(toteEwald,rij,QMD%uv,2.0d0*M_PI*QMD%bv,QMD%omega)
  call calcForceEwald(frcEwald,rij,QMD%uv,2.0d0*M_PI*QMD%bv,QMD%omega)
  call calcStressEwald(strsEwald,rij,QMD%uv,2.0d0*M_PI*QMD%bv,QMD%omega)
  toteEwald=qti*qtj*toteEwald
  frcEwald=-qti*qtj*frcEwald
  strsEwald=qti*qtj*strsEwald
  return
end subroutine Tsuneyuki_Ewald

subroutine Tsuneyuki_Param(z,a,b,c,q,qt)
  implicit none
  real*8 dn
  parameter (dn=0.6d0)
  real*8 QO,aO,bO,cO
  parameter (QO=-1.200d0)
  parameter (aO=2.04740d0)
  parameter (bO=0.17566d0)
  parameter (cO=70.37d0)
  real*8 QSi,aSi,bSi,cSi
  parameter (QSi=2.400d0)
  parameter (aSi=0.8688d0)
  parameter (bSi=0.03285d0)
  parameter (cSi=23.18d0)
!
  real(8):: a,b,c,q,qt
  integer :: z

  if (z==8) then
     a=aO
     b=bO
     c=cO
!     q=QO
     q=-2.d0*dn
     qt=-(1.d0+dn)
  elseif (z==14) then
     a=aSi
     b=bSi
     c=cSi
!     q=QSi
     q=4.d0*dn
     qt=4.d0*dn
  endif   

  return
end subroutine Tsuneyuki_Param

subroutine Tsuneyuki_gSiO(r,eta,gij,dgij)
  implicit none
  real*8 r,eta,gij,dgij
  gij=(1.0d0+eta*r)*exp(-2.0d0*eta*r)
  dgij=-2.0d0*eta*(1.0d0+eta*r)*exp(-2.0d0*eta*r)+eta*exp(-2.0d0*eta*r)
  return
end subroutine Tsuneyuki_gSiO

subroutine Tsuneyuki_gOO(r,eta,gij,dgij)
  implicit none
  real*8 r,eta,gij,dgij
  gij=(1.0d0+11.0d0*(eta*r)/8.0d0+3.0d0*(eta*r)**2/4.0d0+(eta*r)**3/6.0d0) &
     *exp(-2.0d0*eta*r)
  dgij=-2.0d0*eta*(1.0d0+11.0d0*(eta*r)/8.0d0+3.0d0*(eta*r)**2/4.0d0 &
     +(eta*r)**3/6.0d0)*exp(-2.0d0*eta*r) &
     +(11.0d0*eta/8.0d0+3.0d0*eta**2*r/2.0d0+eta**3*r**2/2.0d0)*exp(-2.0d0*eta*r)
  return
end subroutine Tsuneyuki_gOO
!--------1---------2---------3---------4---------5---------6---------7--
!--------1---------2---------3---------4---------5---------6---------7--
subroutine ZRL
! Billeter, Curioni, Fischer and Andreoni, PRB73, 155329 (2006); 
! PRB79, 169904(E) (2009).
! Augmented Tersoff potenntial for Si-O
  use periodic_lattice_module, only: cell_list_type, kr, periodic_lattice_type
  integer :: ipair,i,j,k,l
  integer :: jj
  real*8 :: tote,frc(3,QMD%natom),strs(3,3)
  real*8 :: rri(3),rrj(3)
  real*8 :: rrij(3),raij(3),dij
  type(periodic_lattice_type) :: lattice
  type(cell_list_type), allocatable :: cell_of_replica(:,:)
  real*8 :: dbdri(3,QMD%natom),dbdrj(3,QMD%natom),dbdrk(3,QMD%natom)
  real*8 :: dbdrk_raik(3,3,QMD%natom)
  real*8 :: ddijdri(3),ddijdrj(3),dvijdri(3),dvijdrj(3)
! generalized Morse potential
  real*8 :: totem,frcm(3,QMD%natom),strsm(3,3)
  real*8 :: fijij,dfijijddij,bijij,vr,va,vij
  real*8 :: dvrdd,dvadd,dvadb
! core energies
  real*8 :: totec,frcc(3,QMD%natom),strsc(3,3),e0i,decidzi,decidb
  real*8 :: ddzidd(QMD%natom),ddzidb(QMD%natom),dzidri2(3),dzidrj2(3)
! additional penalty
  real*8 :: totep,frcp(3,QMD%natom),strsp(3,3),eci
  real*8 :: zi,ziz0i,fsz,dfszdz,ci1,ci2,dzi,ddzidzi

  lattice%direct_lattice=real(QMD%uv,kr)
  allocate(cell_of_replica(QMD%natom,QMD%natom))

!$omp parallel do private(j)
  do i=1,QMD%natom
     do j=1,QMD%natom
        cell_of_replica(j,i)=lattice%get_cell_of_replica( &
              real(QMD%rr(:,j),kr), &
              real(QMD%rr(:,i),kr), &
              real(ZRL_sij(QMD%zatm(i),QMD%zatm(j)),kr))
     enddo
  enddo

! generalized Morse potential, eq.(E-13) with eq.(1)
  totem=0.d0
  frcm=0.d0
  strsm=0.d0
  do i=1,QMD%natom
     rri=QMD%rr(:,i)
     do j=1,QMD%natom
        do jj=1,size(cell_of_replica(j,i)%list,dim=2)
           if (j==i.and.all(cell_of_replica(j,i)%list(:,jj)==floor(rri(:)))) cycle
           rrj=cell_of_replica(j,i)%list(:,jj)+modulo(QMD%rr(:,j),1.0d0)
           rrij=rrj(:)-rri(:)
           raij=matmul(QMD%uv,rrij) ! relative to Cartesian
           dij=sqrt(sum(raij**2)) ! bond length in Bohr
           call ZRL_fijij(i,j,dij,fijij,dfijijddij,ipair)
           if (ipair==0) cycle
           call ZRL_bijij(i,j,cell_of_replica(j,i)%list(:,jj),bijij,dbdri,dbdrj,dbdrk,dbdrk_raik)
           call ZRL_VR(i,j,dij,fijij,dfijijddij,vr,dvrdd)
           call ZRL_VA(i,j,dij,fijij,dfijijddij,bijij,va,dvadd,dvadb)
           vij=vr+va
           totem=totem+vij*0.5d0
           call ddijdr(dij,raij,ddijdri,ddijdrj)
           dvijdri=(dvrdd+dvadd)*ddijdri
           dvijdrj=(dvrdd+dvadd)*ddijdrj
           frcm(:,i)=frcm(:,i)-0.5d0*dvijdri
           frcm(:,j)=frcm(:,j)-0.5d0*dvijdrj
           do l=1,3
              strsm(:,l)=strsm(:,l)-0.5d0*dvijdrj*raij(l)
           enddo
           do k=1,QMD%natom
              frcm(:,i)=frcm(:,i)-0.5d0*dvadb*dbdri(:,k)
              frcm(:,j)=frcm(:,j)-0.5d0*dvadb*dbdrj(:,k)
              frcm(:,k)=frcm(:,k)-0.5d0*dvadb*dbdrk(:,k)
              do l=1,3
                 strsm(:,l)=strsm(:,l)-0.5d0*dvadb*dbdrj(:,k)*raij(l) &
                       -0.5d0*dvadb*dbdrk_raik(:,l,k)
              enddo ! l
           enddo ! k 
        enddo ! jj
     enddo ! j
  enddo ! i

! core energies, the second term in eq.(E-13)
  totec=0.d0
  frcc=0.d0
  strsc=0.d0
  do i=1,QMD%natom
     e0i=ZRL_e0i(QMD%zatm(i))
     totec=totec+e0i
  enddo ! i   

! additional penalty, the third term in eq.(E-13) and eq.(14)-(17)
  totep=0.d0
  frcp=0.d0
  strsp=0.d0
  do i=1,QMD%natom
     rri=QMD%rr(:,i)
     zi=0.0d0
     do j=1,QMD%natom
        do jj=1,size(cell_of_replica(j,i)%list,dim=2)
           if (j==i.and.all(cell_of_replica(j,i)%list(:,jj)==floor(rri(:)))) cycle
           rrj=cell_of_replica(j,i)%list(:,jj)+modulo(QMD%rr(:,j),1.0d0)
           rrij=rrj(:)-rri(:)
           raij=matmul(QMD%uv,rrij) ! relative to Cartesian
           dij=sqrt(sum(raij**2)) ! bond length in Bohr
           call ZRL_fijij(i,j,dij,fijij,dfijijddij,ipair)
           if (ipair==0) cycle
           call ZRL_bijij(i,j,cell_of_replica(j,i)%list(:,jj),bijij,dbdri,dbdrj,dbdrk,dbdrk_raik)
           zi=zi+fijij*bijij
        enddo ! jj
     enddo ! j
     ziz0i=zi-ZRL_z0i(QMD%zatm(i))
! eq.(15)(17)
     call ZRL_fsz(i,abs(ziz0i),fsz,dfszdz)
     dzi=ziz0i/abs(ziz0i)*fsz
     ddzidzi=dfszdz

     ci1=ZRL_ci1(QMD%zatm(i))
     ci2=ZRL_ci2(QMD%zatm(i))
! eq.(14)
     eci=ci1*dzi+ci2*dzi**2
     totep=totep+eci
     do j=1,QMD%natom
        do jj=1,size(cell_of_replica(j,i)%list,dim=2)
           if (j==i.and.all(cell_of_replica(j,i)%list(:,jj)==floor(rri(:)))) cycle
           rrj=cell_of_replica(j,i)%list(:,jj)+modulo(QMD%rr(:,j),1.0d0)
           rrij=rrj(:)-rri(:)
           raij=matmul(QMD%uv,rrij) ! relative to Cartesian
           dij=sqrt(sum(raij**2)) ! bond length in Bohr
           call ZRL_fijij(i,j,dij,fijij,dfijijddij,ipair)
           if (ipair==0) cycle
           call ZRL_bijij(i,j,cell_of_replica(j,i)%list(:,jj),bijij,dbdri,dbdrj,dbdrk,dbdrk_raik)

           call ddijdr(dij,raij,ddijdri,ddijdrj)
           dzidri2=bijij*dfijijddij*ddijdri
           dzidrj2=bijij*dfijijddij*ddijdrj

           decidzi=(ci1+2.0d0*ci2*dzi)*ddzidzi
           decidb=decidzi*fijij
           frcp(:,i)=frcp(:,i)-decidzi*dzidri2
           frcp(:,j)=frcp(:,j)-decidzi*dzidrj2
           do l=1,3
              strsp(:,l)=strsp(:,l)-decidzi*dzidrj2*raij(l)
           enddo
           do k=1,QMD%natom
              frcp(:,i)=frcp(:,i)-decidb*dbdri(:,k)
              frcp(:,j)=frcp(:,j)-decidb*dbdrj(:,k)
              frcp(:,k)=frcp(:,k)-decidb*dbdrk(:,k)
              do l=1,3
                 strsp(:,l)=strsp(:,l)-decidb*dbdrj(:,k)*raij(l) &
                                      -decidb*dbdrk_raik(:,l,k)
              enddo ! l
           enddo ! k
        enddo ! jj
     enddo ! j
  enddo ! i   

  deallocate(cell_of_replica)

  QMD%tote=totem+totec+totep
!  QMD%tote=totem+totec
!  QMD%tote=totem
!  QMD%tote=totec
!  QMD%tote=totep
  QMD%frc=frcm+frcc+frcp
!  QMD%frc=frcm+frcc
!  QMD%frc=frcm
!  QMD%frc=frcc
!  QMD%frc=frcp
  QMD%strs=(strsm+strsc+strsp)/QMD%omega

end subroutine ZRL

subroutine ZRL_fijij(i,j,dij,fijij,dfijijddij,ipair)
! eq.(3)(4b) in PRB73,155329(2006)
   implicit none
   real(8) :: dij,fijij,dfijijddij,pi
   real(8) :: ri,rj,rij,si,sj,sij,rtmp
   integer :: ipair,i,j,zi,zj

   pi=4d0*atan(1.0d0)
!   dijang=dij*bohr

! eq.(4b) in PRB73,155329(2006)
   zi=QMD%zatm(i)
   zj=QMD%zatm(j)

   ri=ZRL_ri(zi)
   rj=ZRL_ri(zj)
   rij=sqrt(ri*rj)

   si=ZRL_si(zi)
   sj=ZRL_si(zj)
   sij=sqrt(si*sj)

! eq.(3) in PRB73,155329(2006)
   if (dij<=rij) then
      fijij=1.0d0
      dfijijddij=0.0d0
      ipair=1
   elseif (dij<=sij) then
!      rtmp=(dij-rij)/(sij-rij)
      rtmp=pi*(dij-rij)/(sij-rij)
      fijij=0.5d0*(1.0d0+cos(rtmp))
      dfijijddij=-0.5*sin(rtmp)*pi/(sij-rij)
      ipair=1
   else
      fijij=0.d0
      dfijijddij=0.0d0
      ipair=0
   endif   
end subroutine ZRL_fijij

subroutine ZRL_bijij(i,j,shiftj,bijij,dbdri,dbdrj,dbdrk,dbdrk_raik)
! eq.(6) in PRB73,155329(2006)
   use periodic_lattice_module, only: cell_list_type, kr, periodic_lattice_type
   implicit none
   real(8) :: bijij,pi
   real(8) :: dbdri(3,QMD%natom),dbdrj(3,QMD%natom),dbdrk(3,QMD%natom)
   real(8) :: dbdrk_raik(3,3,QMD%natom)
   real(8) :: rri(3),rrj(3),rrk(3),rrij(3),rrik(3),raij(3),raik(3)
   integer :: shiftj(3)
   type(periodic_lattice_type) :: lattice
   type(cell_list_type) :: cell_of_replica
   real(8) :: dij,dik,cosjik,fikik,eijkijk,tijki,ztijij,chiij,betai,ni
   real(8) :: dfikikddik,deddij,deddik,dbijijdztij
   real(8) :: dcdri(3),dcdrj(3),dcdrk(3)
   real(8) :: ddijdri(3),ddijdrj(3),ddikdri(3),ddikdrk(3)
   real(8) :: dztijijddij,dztijijddik,dztijijdcos,dtijkidcos
   real(8) :: tempk(3)
   integer :: ipair,i,j,k,zi,zj
   integer :: kk,l

!   pi=4.0d0*atan(1.0d0)
   zi=QMD%zatm(i)
   zj=QMD%zatm(j)
   chiij=ZRL_chi(zi,zj)
   betai=ZRL_betai(zi)
   ni=ZRL_ni(zi)

   dbdri=0.0d0
   dbdrj=0.0d0
   dbdrk=0.0d0
   ztijij=0.0d0
   dbdrk_raik=0.0d0
   rri=QMD%rr(:,i)
   rrj=shiftj(:)+modulo(QMD%rr(:,j),1.0d0)
   rrij=rrj(:)-rri(:)
   raij=matmul(QMD%uv,rrij) ! relative to Cartesian
   dij=sqrt(sum(raij**2)) ! bond length in Bohr

   lattice%direct_lattice=real(QMD%uv,kr)

   do k=1,QMD%natom
      cell_of_replica=lattice%get_cell_of_replica(real(QMD%rr(:,k),kr),real(QMD%rr(:,i),kr),real(ZRL_sij(zi,QMD%zatm(k)),kr))
      do kk=1,size(cell_of_replica%list,dim=2)
         if (k==i.and.all(cell_of_replica%list(:,kk)==floor(rri(:)))) cycle
         if (k==j.and.all(cell_of_replica%list(:,kk)==shiftj(:))) cycle
         rrk=cell_of_replica%list(:,kk)+modulo(QMD%rr(:,k),1.0d0)
         rrik=rrk(:)-rri(:)
         raik=matmul(QMD%uv,rrik) ! relative to Cartesian
         dik=sqrt(sum(raik**2)) ! bond length in Bohr

         cosjik=dot_product(raij,raik)/(dij*dik)

         call ZRL_fijij(i,k,dik,fikik,dfikikddik,ipair) ! eq.(3)
         if (ipair==0) cycle
         call ZRL_eijkijk(i,j,k,dij,dik,eijkijk,deddij,deddik) ! eq.(11)
         call ZRL_tijki(i,cosjik,tijki,dtijkidcos) ! eq.(8)
         ztijij=ztijij+fikik*eijkijk*tijki ! eq.(7)
         dztijijddij=fikik*deddij*tijki
         dztijijddik=fikik*deddik*tijki+dfikikddik*eijkijk*tijki
         dztijijdcos=fikik*eijkijk*dtijkidcos

         call ddijdr(dij,raij,ddijdri,ddijdrj)
         call ddijdr(dik,raik,ddikdri,ddikdrk)
         call dcosjikdr(raij,raik,dij,dik,cosjik,dcdri,dcdrj,dcdrk)
         dbdri(:,k)=dbdri(:,k)+dztijijddij*ddijdri+dztijijddik*ddikdri+dztijijdcos*dcdri
         dbdrj(:,k)=dbdrj(:,k)+dztijijddij*ddijdrj                    +dztijijdcos*dcdrj
         tempk     =                              +dztijijddik*ddikdrk+dztijijdcos*dcdrk
         dbdrk(:,k)=dbdrk(:,k)+tempk(:)

         do l=1,3
            dbdrk_raik(:,l,k)=dbdrk_raik(:,l,k)+tempk(:)*raik(l)
         enddo ! l
      enddo ! kk
   enddo ! k   

   bijij=chiij*(1.0d0+(betai*ztijij)**ni)**(-0.5d0/ni)
!   bijij=chiij*(1.0d0+(betai*ztijij)**ni)**(-0.5d0*ni)
   if (ztijij==0.0d0) then
      dbdri=0.0d0
      dbdrj=0.0d0
      dbdrk=0.0d0
      dbdrk_raik=0.0d0
   else
      dbijijdztij=-0.5d0/ni*chiij*(1.0d0+(betai*ztijij)**ni)**(-0.5d0/ni-1.0d0)*ni*(betai*ztijij)**(ni-1.0d0)*betai
      dbdri=dbijijdztij*dbdri
      dbdrj=dbijijdztij*dbdrj
      dbdrk=dbijijdztij*dbdrk
      dbdrk_raik=dbijijdztij*dbdrk_raik
   endif

end subroutine ZRL_bijij

subroutine ZRL_eijkijk(i,j,k,dij,dik,eijkijk,deddij,deddik)
! eq.(11) in PRB73,155329(2006)
   implicit none
   real(8) :: dij,dik,eijkijk,deddij,deddik,mui,muj,muk,muij,muik
   integer:: i,j,k,zi,zj,zk,mi

   zi=QMD%zatm(i)
   zj=QMD%zatm(j)
   zk=QMD%zatm(k)
! eq.(5) and Table3 in PRB73,155329(2006)
! Table2 in PRB79,169904(E)(2009)
   muij=ZRL_muij(zi,zj)
   muik=ZRL_muij(zi,zk)
! Table2 in PRB79,169904(E)(2009)
   mi=ZRL_mi(zi)
! eq.(11) in PRB73,155329(2006)
   eijkijk=exp((muij*dij-muik*dik)**mi)
   deddij= eijkijk*(muij*dij-muik*dik)**(mi-1)*mi*muij
   deddik=-eijkijk*(muij*dij-muik*dik)**(mi-1)*mi*muik
end subroutine ZRL_eijkijk

subroutine ZRL_tijki(i,cosjik,tijki,dtijkidcos)
! eq.(8) in PRB73,155329(2006)
   implicit none
   real(8) :: cosjik,tijki,dtijkidcos,ci,di,hi
   integer:: i,zi
! Table2 in PRB79,169904(E)(2009)
   zi=QMD%zatm(i)
   ci=ZRL_ci(zi)
   di=ZRL_di(zi)
   hi=ZRL_hi(zi)
! eq.(8) in PRB73,155329(2006)
   tijki=1.0d0+(ci/di)**2-ci**2/(di**2+(hi-cosjik)**2)
   dtijkidcos=-ci**2/(di**2+(hi-cosjik)**2)**2*2.0d0*(hi-cosjik)
end subroutine ZRL_tijki

subroutine ZRL_VR(i,j,dij,fijij,dfijijddij,vr,dvrdd)
! eq.(8) in PRB73,155329(2006)
   implicit none
   real(8) :: dij,fijij,dfijijddij,vr,dvrdd,ai,aj,aij,lmdi,lmdj,lmdij
   integer:: i,j,zi,zj
   zi=QMD%zatm(i)
   zj=QMD%zatm(j)
! eq.(4)(5) and Table3 in PRB73,155329(2006)
! Table2 in PRB79,169904(E)(2009)
   aij=ZRL_aij(zi,zj)
   lmdij=ZRL_lmdij(zi,zj)
! eq.(2) in PRB73,155329(2006)
   vr=fijij*aij*exp(-lmdij*dij)
   dvrdd=vr*(-lmdij)+dfijijddij*aij*exp(-lmdij*dij)
end subroutine ZRL_VR

subroutine ZRL_VA(i,j,dij,fijij,dfijijddij,bijij,va,dvadd,dvadb)
! eq.(8) in PRB73,155329(2006)
   implicit none
   real(8) :: dij,fijij,dfijijddij,bijij,va,dvadd,dvadb,bi,bj,bij,mui,muj,muij
   integer:: i,j,zi,zj
   zi=QMD%zatm(i)
   zj=QMD%zatm(j)
! eq.(4)(5) and Table3 in PRB73,155329(2006)
! Table2 in PRB79,169904(E)(2009)
   bij=ZRL_bij(zi,zj)
   muij=ZRL_muij(zi,zj)
! eq.(2) in PRB73,155329(2006)
   va=-fijij*bijij*bij*exp(-muij*dij)
   dvadd=va*(-muij)-dfijijddij*bijij*bij*exp(-muij*dij)
   dvadb=-fijij*bij*exp(-muij*dij)
end subroutine ZRL_VA

subroutine ZRL_fsz(i,z,fsz,dfszdz)
! eq.(17) in PRB73,155329(2006)
   implicit none
   real(8) :: z,fsz,dfszdz,zt,zb,pi,rtmp,x
   integer:: i
   pi=4.0d0*atan(1.0d0)
! Table4
   zt=0.49751
   zb=0.20039
!   az=abs(z)
! eq.(17)
   if (z-floor(z)<=zt-zb) then
      fsz=0.0d0
      dfszdz=0.0d0
   elseif (z-floor(z)<=zt+zb) then
      x=0.5d0*(z-floor(z)-zt)/zb
!      rtmp=pi*(x-nint(x))
      rtmp=pi*x
      fsz=0.5d0*(1.0d0+sin(rtmp))
      dfszdz=0.5d0*cos(rtmp)*0.5d0*pi/zb
!      fsz=0.5d0*(1.0d0+sin((z-zt)/zb))
!      dfszdz=0.5d0*cos((z-zt)/zb)*pi/zb
   else
      fsz=1.0d0
      dfszdz=0.0d0
   endif   
   fsz=fsz+floor(z)
end subroutine ZRL_fsz

subroutine ZRL_fsz2(i,z,fsz,dfszdz)
! eq.(17) in PRB73,155329(2006)
   implicit none
   real(8) :: z,fsz,dfszdz,zt,zb,pi,rtmp,x
   integer:: i
   pi=4.0d0*atan(1.0d0)
! Table4
   zt=0.49751
   zb=0.20039
!   az=abs(z)
! eq.(17)
!   if (z<=zt-zb) then
!      fsz=0.0d0
!      dfszdz=0.0d0
!   elseif (z<=zt+zb) then
      x=(z-zt)/zb
      rtmp=pi*(x-nint(x))
      fsz=0.5d0*(1.0d0+sin(rtmp))+nint(x)
      dfszdz=0.5d0*cos(rtmp)*pi/zb
!!      fsz=0.5d0*(1.0d0+sin((z-zt)/zb))
!!      dfszdz=0.5d0*cos((z-zt)/zb)*pi/zb
!   else
!      fsz=1.0d0
!      dfszdz=0.0d0
!   endif   
!   fsz=fsz+int(z)
end subroutine ZRL_fsz2

real(8) function ZRL_ai(z)
  implicit none
  integer :: z
  if (z==14) then
     ZRL_ai=1830.79
  elseif (z==8) then
     ZRL_ai=3331.06
  elseif (z==7) then
     ZRL_ai=6368.21
  elseif (z==1) then
     ZRL_ai=86.9235
  else
     ZRL_ai=0.d0
  endif
  ZRL_ai=ZRL_ai/hartree ! eV to Hartree
end function ZRL_ai

real(8) function ZRL_aij(zi,zj)
  implicit none
  integer :: zi,zj
  real(8)::ai,aj
  ai=ZRL_ai(zi)
  aj=ZRL_ai(zj)
  ZRL_aij=sqrt(ai*aj)
  if ((zi==14.and.zj==8).or.(zj==14.and.zi==8)) then
     ZRL_aij=ZRL_aij*1.04752
  elseif ((zi==14.and.zj==7).or.(zj==14.and.zi==7)) then   
     ZRL_aij=ZRL_aij*0.58647
  elseif ((zi==14.and.zj==1).or.(zj==14.and.zi==1)) then   
     ZRL_aij=ZRL_aij*1.52966
  elseif ((zi==8.and.zj==7).or.(zj==8.and.zi==7)) then   
     ZRL_aij=ZRL_aij*1.26527
  elseif ((zi==8.and.zj==1).or.(zj==8.and.zi==1)) then   
     ZRL_aij=ZRL_aij*0.99853
  elseif ((zi==7.and.zj==1).or.(zj==7.and.zi==1)) then   
     ZRL_aij=ZRL_aij*0.83424
  endif
end function ZRL_aij

real(8) function ZRL_bi(z)
  implicit none
  integer :: z
  if (z==14) then
     ZRL_bi=471.195
  elseif (z==8) then
     ZRL_bi=260.477
  elseif (z==7) then
     ZRL_bi=511.205
  elseif (z==1) then
     ZRL_bi=42.9815
  else
     ZRL_bi=0.d0
  endif
  ZRL_bi=ZRL_bi/hartree ! eV to Hartree
end function ZRL_bi

real(8) function ZRL_bij(zi,zj)
  implicit none
  integer :: zi,zj
  real(8)::bi,bj
  bi=ZRL_bi(zi)
  bj=ZRL_bi(zj)
  ZRL_bij=sqrt(bi*bj)
  if ((zi==14.and.zj==8).or.(zj==14.and.zi==8)) then
     ZRL_bij=ZRL_bij*0.99978
  elseif ((zi==14.and.zj==7).or.(zj==14.and.zi==7)) then   
     ZRL_bij=ZRL_bij*1.10293
  elseif ((zi==14.and.zj==1).or.(zj==14.and.zi==1)) then   
     ZRL_bij=ZRL_bij*1.68173
  elseif ((zi==8.and.zj==7).or.(zj==8.and.zi==7)) then   
     ZRL_bij=ZRL_bij*1.00075
  elseif ((zi==8.and.zj==1).or.(zj==8.and.zi==1)) then   
     ZRL_bij=ZRL_bij*1.01274
  elseif ((zi==7.and.zj==1).or.(zj==7.and.zi==1)) then   
     ZRL_bij=ZRL_bij*0.97237
  endif
end function ZRL_bij

real(8) function ZRL_lmdi(z)
  implicit none
  integer :: z
  if (z==14) then
     ZRL_lmdi=2.62392
  elseif (z==8) then
     ZRL_lmdi=3.78336
  elseif (z==7) then
     ZRL_lmdi=5.60181
  elseif (z==1) then
     ZRL_lmdi=3.8593
  else
     ZRL_lmdi=0.d0
  endif
  ZRL_lmdi=ZRL_lmdi*bohr ! Ang^{-1} to Bohr^{-1}
end function ZRL_lmdi

real(8) function ZRL_lmdij(zi,zj)
  implicit none
  integer :: zi,zj
  real(8)::lmdi,lmdj
  lmdi=ZRL_lmdi(zi)
  lmdj=ZRL_lmdi(zj)
  ZRL_lmdij=(lmdi+lmdj)*0.5d0
  if ((zi==14.and.zj==8).or.(zj==14.and.zi==8)) then
     ZRL_lmdij=ZRL_lmdij+0.74003*bohr
  elseif ((zi==14.and.zj==7).or.(zj==14.and.zi==7)) then   
     ZRL_lmdij=ZRL_lmdij-0.73787*bohr
  elseif ((zi==14.and.zj==1).or.(zj==14.and.zi==1)) then   
     ZRL_lmdij=ZRL_lmdij-0.15903*bohr
  elseif ((zi==8.and.zj==7).or.(zj==8.and.zi==7)) then   
     ZRL_lmdij=ZRL_lmdij+2.34383*bohr
  elseif ((zi==8.and.zj==1).or.(zj==8.and.zi==1)) then   
     ZRL_lmdij=ZRL_lmdij+1.03160*bohr
  elseif ((zi==7.and.zj==1).or.(zj==7.and.zi==1)) then   
     ZRL_lmdij=ZRL_lmdij+0.07480*bohr
  endif
end function ZRL_lmdij

real(8) function ZRL_mui(z)
  implicit none
  integer :: z
  if (z==14) then
     ZRL_mui=1.88891
  elseif (z==8) then
     ZRL_mui=3.34402
  elseif (z==7) then
     ZRL_mui=3.16170
  elseif (z==1) then
     ZRL_mui=1.97047
  else
     ZRL_mui=0.d0
  endif
  ZRL_mui=ZRL_mui*bohr ! Ang^{-1} to Bohr^{-1}
end function ZRL_mui

real(8) function ZRL_muij(zi,zj)
  implicit none
  integer :: zi,zj
  real(8)::mui,muj
  mui=ZRL_mui(zi)
  muj=ZRL_mui(zj)
  ZRL_muij=(mui+muj)*0.5d0
  if ((zi==14.and.zj==8).or.(zj==14.and.zi==8)) then
     ZRL_muij=ZRL_muij-0.30051*bohr
  elseif ((zi==14.and.zj==7).or.(zj==14.and.zi==7)) then   
     ZRL_muij=ZRL_muij-0.19843*bohr
  elseif ((zi==14.and.zj==1).or.(zj==14.and.zi==1)) then   
     ZRL_muij=ZRL_muij+0.22168*bohr
  elseif ((zi==8.and.zj==7).or.(zj==8.and.zi==7)) then   
     ZRL_muij=ZRL_muij+3.50573*bohr
  elseif ((zi==8.and.zj==1).or.(zj==8.and.zi==1)) then   
     ZRL_muij=ZRL_muij-0.22005*bohr
  elseif ((zi==7.and.zj==1).or.(zj==7.and.zi==1)) then   
     ZRL_muij=ZRL_muij+0.21563*bohr
  endif
end function ZRL_muij

real(8) function ZRL_ri(z)
  implicit none
  integer :: z
  if (z==14) then
     ZRL_ri=2.44809
  elseif (z==8) then
     ZRL_ri=2.26069
  elseif (z==7) then
     ZRL_ri=1.75256
  elseif (z==1) then
     ZRL_ri=0.77985
  else
     ZRL_ri=0.d0
  endif
  ZRL_ri=ZRL_ri/bohr ! Ang to Bohr
end function ZRL_ri

real(8) function ZRL_si(z)
  implicit none
  integer :: z
  if (z==14) then
     ZRL_si=3.08354
  elseif (z==8) then
     ZRL_si=3.31294
  elseif (z==7) then
     ZRL_si=2.41523
  elseif (z==1) then
     ZRL_si=0.88641
  else
     ZRL_si=0.d0
  endif
  ZRL_si=ZRL_si/bohr ! Ang to Bohr
end function ZRL_si

real(8) function ZRL_sij(zi,zj)
   implicit none
   integer, intent(in) :: zi,zj
   ZRL_sij=sqrt(ZRL_si(zi)*ZRL_si(zj))
end function ZRL_sij

real(8) function ZRL_betai(z)
  implicit none
  integer :: z
  if (z==14) then
     ZRL_betai=1.0999d-6
  elseif (z==8) then
     ZRL_betai=1.0027
  elseif (z==7) then
     ZRL_betai=4.4422d-3
  elseif (z==1) then
     ZRL_betai=4.0d0
  else
     ZRL_betai=0.d0
  endif
end function ZRL_betai

real(8) function ZRL_ni(z)
  implicit none
  integer :: z
  if (z==14) then
     ZRL_ni=0.78766
  elseif (z==8) then
     ZRL_ni=3.98638
  elseif (z==7) then
     ZRL_ni=2.42635
  elseif (z==1) then
     ZRL_ni=1.00921
  else
     ZRL_ni=0.d0
  endif
end function ZRL_ni

integer function ZRL_mi(z)
  implicit none
  integer :: z
  if (z==14) then
     ZRL_mi=3
  elseif (z==8) then
     ZRL_mi=1
  elseif (z==7) then
     ZRL_mi=1
  elseif (z==1) then
     ZRL_mi=1
  else
     ZRL_mi=0
  endif
end function ZRL_mi

real(8) function ZRL_ci(z)
  implicit none
  integer :: z
  if (z==14) then
     ZRL_ci=1.0039d5
  elseif (z==8) then
     ZRL_ci=0.0d0
  elseif (z==7) then
     ZRL_ci=2.2955d4
  elseif (z==1) then
     ZRL_ci=0.0d0
  else
     ZRL_ci=0.d0
  endif
end function ZRL_ci

real(8) function ZRL_di(z)
  implicit none
  integer :: z
  if (z==14) then
     ZRL_di=16.21701
  elseif (z==8) then
     ZRL_di=1.0d0
  elseif (z==7) then
     ZRL_di=24.78674
  elseif (z==1) then
     ZRL_di=1.0d0
  else
     ZRL_di=0.d0
  endif
end function ZRL_di

real(8) function ZRL_hi(z)
  implicit none
  integer :: z
  if (z==14) then
     ZRL_hi=-0.59784
  elseif (z==8) then
     ZRL_hi=-0.52909
  elseif (z==7) then
     ZRL_hi=-0.45450
  elseif (z==1) then
     ZRL_hi=0.96783
  else
     ZRL_hi=0.d0
  endif
end function ZRL_hi

real(8) function ZRL_chi(zi,zj)
  implicit none
  integer :: zi,zj
  ZRL_chi=1.0d0
end function ZRL_chi

real(8) function ZRL_e0i(z)
  implicit none
  integer :: z
  if (z==14) then
     ZRL_e0i=-103.733
  elseif (z==8) then
     ZRL_e0i=-432.158
  elseif (z==7) then
     ZRL_e0i=-264.156
  elseif (z==1) then
     ZRL_e0i=-13.174
  else
     ZRL_e0i=0.d0
  endif
  ZRL_e0i=ZRL_e0i/hartree ! eV to Hartree
end function ZRL_e0i

real(8) function ZRL_z0i(z)
  implicit none
  integer :: z
  if (z==14) then
     ZRL_z0i=3.70
  elseif (z==8) then
     ZRL_z0i=2.80
  elseif (z==7) then
     ZRL_z0i=1.75
  elseif (z==1) then
     ZRL_z0i=1.00
  else
     ZRL_z0i=0.d0
  endif
end function ZRL_z0i

real(8) function ZRL_ci1(z)
  implicit none
  integer :: z
  if (z==14) then
     ZRL_ci1=-0.1238
  elseif (z==8) then
     ZRL_ci1=-0.0038
  elseif (z==7) then
     ZRL_ci1=-0.0868
  elseif (z==1) then
     ZRL_ci1=0.0d0
  else
     ZRL_ci1=0.d0
  endif
  ZRL_ci1=ZRL_ci1/hartree ! eV to Hartree
end function ZRL_ci1

real(8) function ZRL_ci2(z)
  implicit none
  integer :: z
  if (z==14) then
     ZRL_ci2=0.2852
  elseif (z==8) then
     ZRL_ci2=0.1393
  elseif (z==7) then
     ZRL_ci2=0.2454
  elseif (z==1) then
     ZRL_ci2=0.0d0
  else
     ZRL_ci2=0.d0
  endif
  ZRL_ci2=ZRL_ci2/hartree ! eV to Hartree
end function ZRL_ci2

subroutine ADP_KWU14
  use adp_kwu14_module, only: adp_kwu14_parameter, adp_kwu14_type
  implicit none

  real(8) :: uv(3,3)
  type(adp_kwu14_type) :: pot

  uv=QMD%uv*bohr ! Bohr to Ang

  pot%parameter = adp_kwu14_parameter
  call pot%set_system(uv,QMD%zatm,QMD%rr)

  QMD%tote=pot%energy() ! in eV
  QMD%tote=QMD%tote/hartree ! eV to Hartree

  QMD%frc=pot%force() ! in eV/Ang
  QMD%frc=QMD%frc/hartree*bohr ! eV/Ang to Hartree/Bohr

  QMD%strs=pot%stress() ! in eV/Ang^3
  QMD%strs=QMD%strs/hartree*bohr**3 ! eV/Ang^3 to Hartree/Bohr^3

end subroutine ADP_KWU14
!--------1---------2---------3---------4---------5---------6---------7--
!--------1---------2---------3---------4---------5---------6---------7--
subroutine Jmatgen
  use m_jmatgen, only: calcJmatgen

  integer :: natom
  real*8 :: cell(3,3)
  character(len=4) :: atom(QMD%natom)
  real*8 :: coord(3,QMD%natom)
  real*8 :: energy
  real*8 :: force(3,QMD%natom)
  real*8 :: virial(3,3)
  real*8 :: stress(3,3)

  integer :: ia

! initialize
  natom = QMD%natom
  cell = QMD%uv*bohr ! Bohr to Ang
  do ia=1,natom
    atom = getAtomicName(QMD%zatm(ia)) ! Element number to Element name
  enddo
  coord = QMD%ra*bohr  ! Bohr to Ang
  energy = 0.d0 ! eV
  force = 0.d0 ! Hartree/Bohr
  virial = 0.d0 ! Hartree

! calculate
  call calcJmatgen(natom, cell, atom, coord, energy, force, virial)
  stress = virial/QMD%omega ! Hartree/Bohr^3

! total
  QMD%tote = energy/hartree ! eV to Hartree
  QMD%frc = force ! Hartree/Bohr
  QMD%strs = stress ! Hartree/Bohr^3

contains
  !!----------------
  !! functions for element name.
  !!----------------
  function getAtomicName( number ) result(symbol)
    integer, intent(in) :: number
    character(len=4) :: symbol

    select case( number )
    case( 1); symbol= "H"; case( 2); symbol="He"; case( 3); symbol="Li";
    case( 4); symbol="Be"; case( 5); symbol= "B"; case( 6); symbol= "C";
    case( 7); symbol= "N"; case( 8); symbol= "O"; case( 9); symbol= "F";
    case(10); symbol="Ne"; case(11); symbol="Na"; case(12); symbol="Mg";
    case(13); symbol="Al"; case(14); symbol="Si"; case(15); symbol= "P";
    case(16); symbol= "S"; case(17); symbol="Cl"; case(18); symbol="Ar";
    case(19); symbol= "K"; case(20); symbol="Ca"; case(21); symbol="Sc";
    case(22); symbol="Ti"; case(23); symbol= "V"; case(24); symbol="Cr";
    case(25); symbol="Mn"; case(26); symbol="Fe"; case(27); symbol="Co";
    case(28); symbol="Ni"; case(29); symbol="Cu"; case(30); symbol="Zn";
    case(31); symbol="Ga"; case(32); symbol="Ge"; case(33); symbol="As";
    case(34); symbol="Se"; case(35); symbol="Br"; case(36); symbol="Kr";
    case(37); symbol="Rb"; case(38); symbol="Sr"; case(39); symbol= "Y";
    case(40); symbol="Zr"; case(41); symbol="Nb"; case(42); symbol="Mo";
    case(43); symbol="Tc"; case(44); symbol="Ru"; case(45); symbol="Rh";
    case(46); symbol="Pd"; case(47); symbol="Ag"; case(48); symbol="Cd";
    case(49); symbol="In"; case(50); symbol="Sn"; case(51); symbol="Sb";
    case(52); symbol="Te"; case(53); symbol= "I"; case(54); symbol="Xe";
    case(55); symbol="Cs"; case(56); symbol="Ba"; case(57); symbol="La";
    case(58); symbol="Ce"; case(59); symbol="Pr"; case(60); symbol="Nd";
    case(61); symbol="Pm"; case(62); symbol="Sm"; case(63); symbol="Eu";
    case(64); symbol="Gd"; case(65); symbol="Tb"; case(66); symbol="Dy";
    case(67); symbol="Ho"; case(68); symbol="Er"; case(69); symbol="Tm";
    case(70); symbol="Yb"; case(71); symbol="Lu"; case(72); symbol="Hf";
    case(73); symbol="Ta"; case(74); symbol= "W"; case(75); symbol="Re";
    case(76); symbol="Os"; case(77); symbol="Ir"; case(78); symbol="Pt";
    case(79); symbol="Au"; case(80); symbol="Hg"; case(81); symbol="Tl";
    case(82); symbol="Pb"; case(83); symbol="Bi"; case(84); symbol="Po";
    case(85); symbol="At"; case(86); symbol="Rn"; case(87); symbol="Fr";
    case(88); symbol="Ra"; case(89); symbol="Ac"; case(90); symbol="Th";
    case(91); symbol="Pa"; case(92); symbol= "U"; case(93); symbol="Np";
    case(94); symbol="Pu"; case(95); symbol="Am"; case(96); symbol="Cm";
    case(97); symbol="Bk"; case(98); symbol="Cf"; case(99); symbol="Es";
    case(100); symbol="Fm"; case(101); symbol="Md"; case(102); symbol="No";
    case(103); symbol="Lr";
    case default
       write(*,'(a,i8)') '# Error!: unknown element number:', number
       stop
    end select
  end function getAtomicName

end subroutine Jmatgen

end module
