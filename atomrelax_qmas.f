! from atomrelax.f in QMAS
! change the name of variables
! 150610, TM 
!
! Last updated: 20101222
! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**
!  2008.02.07
      subroutine atomrelax
! *                                                                    *
! *     imd=1: simple molecular dynamics                               *
! *        =2: simple ralaxation (Steepest decent scheme)              *
! *        =3: GDIIS                                                   *
! *        =4: quenched MD                                             *
! *        =5: damped MD                                               *
! *        =6: quickmin                                                *
! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**
      use paramlist
      implicit none

      integer mlin,n,isp,il,ina
      real*8 rdmax,scosth,gnorm1m,gnorm0,gnorm,dtoter
      real*8, allocatable:: rrsave(:,:)

      allocate (rrsave(3,natom))
      rrsave=rr

! --- parameters
      mlin=nloopa
      lgdiis=2
      rdmax=0.25d0

! --- setup
      if ((loopa==1).and.(loopc==1)) then
        allocate (rah(3,natom,mlin),gradh(3,natom,mlin))
        allocate (hessi2(natom*3,natom*3))
      end if

      if (loopa==1) then
        igdiis=0
      end if

! --- imod: current step,  imod0: previous step,  imod2: next step
      imod=mod((loopa-1),mlin)+1
      if (loopa>1) then
        imod0=mod((loopa-2),mlin)+1
      else
        imod0=mlin
      end if
      imod2=mod(loopa,mlin)+1

! --- input or output of tstep (file: ***tstep.txt)
      if(loopa==1) then
        call outtstep
      else
        call intstep
      end if

! --- relative --> absolute
      ra=matmul(uv,rr)

! --- Saving of Positions and Graidents at the previous relaxation step
      do ina=1,natom
        do il=1,3
        rah(il,ina,imod)=ra(il,ina)
        gradh(il,ina,imod)=-frc(il,ina)
        end do
      end do


! --- information of the previous relaxation step
      call pre_atomrelax(mlin,scosth,gnorm1m,gnorm0,gnorm,dtoter)

! --- atomic relaxation
      if(iabs(imd).eq.1) then
        call simple_md
      else if(iabs(imd).eq.2) then
        if(loopa.le.3) then
          call quickmin(rdmax,scosth)
        else
          call gdiis(mlin,rdmax,scosth,gnorm1m,gnorm0,gnorm,dtoter)
        end if
      else if(iabs(imd).eq.3) then
        call simple_relax
      else if(iabs(imd).eq.4) then
        call quench_md(rdmax)
      else if(iabs(imd).eq.5) then
        if(loopa.le.3) then
          call quickmin(rdmax,scosth)
        else
          call damped_md(rdmax,scosth)
        end if
      else if(iabs(imd).eq.6) then
        call quickmin(rdmax,scosth)
      else
        call simple_relax
      end if

! --- absolute --> relative
      do ina=1,natom
         rr(:,ina)=matmul(ra(:,ina),bv)
      enddo

      toter0=tote

      do ina=1,natom
        if (iposfix(ina).eq.0) then
          rr(:,ina)=rrsave(:,ina)
        end if
      end do
      ra=matmul(uv,rr)

      deallocate (rrsave)

      return
      end

! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**
      subroutine pre_atomrelax(mlin,scosth,gnorm1m,gnorm0,gnorm,dtoter)
! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**
      use paramlist
      implicit none

      integer i,ina,ika,mlin,imod1m
      real*8 scosth,prd,famp,vamp,gnorm,gnorm0,prdt,dtmp,dtoter,
     & gnorm1m,rfac

      do ina=1,natom
        if(iposfix(ina).eq.1) then
!          ika=katm(ina)
          rfac=tstep*omegai/(mass(ina)*mfac(ina))
          vrr(:,ina)=vrr(:,ina)-rfac*matmul(frc(:,ina),bv)
        endif
      enddo  

! ***** velocity*force for the whole system
      famp=0.0d0
      vamp=0.0d0
      prd=0.0d0
      do ina=1,natom
        if(iposfix(ina).eq.1) then
 !         ika=katm(ina)
          famp=famp+sum(frc(:,ina)**2)
          vamp=vamp+sum(matmul(uv,vrr(:,ina))**2)
        end if
      end do
      famp=dsqrt(famp)
      vamp=dsqrt(vamp)
      if(famp.lt.1.0d-12.or.vamp.lt.1.0d-12) then
        scosth=1.0d0
      else
        scosth=prd/(famp*vamp)
      end if
      write(6,'(1x/1x,"costh for whole system =",d16.8)') scosth

      do ina=1,natom
        if(iposfix(ina).eq.1) then
!          ika=katm(ina)
        rfac=tstep*omegai/(mass(ina)*mfac(ina))
        vrr(:,ina)=vrr(:,ina)+rfac*matmul(frc(:,ina),bv)
        end if
      end do

! ***** norm of gradient in previous steps
      imod1m=mod(loopa-3,mlin)+1

      if(loopa.ge.3) then
        gnorm=0.0d0
        gnorm0=0.0d0
        gnorm1m=0.0d0
        do ina=1,natom
          if(iposfix(ina).eq.1) then
          gnorm=gnorm+sum(gradh(:,ina,imod)**2)
          gnorm0=gnorm0+sum(gradh(:,ina,imod0)**2)
          gnorm1m=gnorm1m+sum(gradh(:,ina,imod1m)**2)
          end if
        end do
        gnorm=dsqrt(gnorm)
        gnorm0=dsqrt(gnorm0)
        gnorm1m=dsqrt(gnorm1m)

        dtoter=toter0-tote
        write(6,'(1x,"gnorm-1 gnorm0 gnorm1 dtoter =",4d12.4)') 
     &    gnorm1m,gnorm0,gnorm,dtoter
      end if

      return
      end

! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**
      subroutine simple_md
! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**

      use paramlist
      implicit none

      integer ina
      real*8 dxyz(3)

      write(6,'(1x/1x,"simple molecular dynamics")')
      write(6,'(1x,"atom    dx          dy          dz"/
     & 1x,"------------------------------------------")')

      do ina=1,natom
        dxyz=tstep*matmul(uv,vrr(:,ina))
        ra(:,ina)=ra(:,ina)+dxyz(:)
        write(6,'(2x,i3,2x,3d12.4)') ina,dxyz(1:3)
      end do

      write(6,'(1x,"------------------------------------------"/)')

      return
      end

! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**
      subroutine simple_relax
! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**

      use paramlist
      implicit none

      integer ina
      real*8 dxyz(3)

!      write(6,'(1x/1x,"simple relaxation")')
!      write(6,'(1x,"atom    dx          dy          dz"/
!     & 1x,"------------------------------------------")')

      do ina=1,natom
        if(iposfix(ina).eq.1) then
        dxyz=tstep*matmul(uv,vrr(:,ina))
        ra(:,ina)=ra(:,ina)+dxyz(:)
!        write(6,'(2x,i3,2x,3d12.4)') ina,dxyz(1:3)
        end if
        vrr(:,ina)=0.0d0
      end do

!      write(6,'(1x,"------------------------------------------"/)')

      return
      end

! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**
      subroutine quench_md(rdmax)
! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**
      use paramlist
      implicit none

      integer i,ina
      real*8 dxyz(3),dd,rdmax,prd

      write(6,'(1x/1x,"queched MD")')
      write(6,'(1x,"atom    dx          dy          dz"/
     & 1x,"------------------------------------------")')

      do ina=1,natom
        if(iposfix(ina).eq.1) then
          prd=sum(matmul(uv,vrr(:,ina))*frc(:,ina))
          if (prd.lt.(-tol)) then
             rah(:,ina,imod2)=rah(:,ina,imod)
             vrr(:,ina)=0.0d0
          end if

          dxyz(:)=tstep*matmul(uv,vrr(:,ina))
          dd=dsqrt(sum(dxyz**2))
          if(dd.gt.rdmax) then
            dxyz=dxyz*rdmax/dd
            dd=rdmax
          end if

          ra(:,ina)=ra(:,ina)+dxyz(:)
          write(6,'(2x,i3,2x,3d12.4)') ina,dxyz(1:3)
        end if
      end do

      write(6,'(1x,"------------------------------------------"/)')

      return
      end

! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**
      subroutine quickmin(rdmax,scosth)
! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**
      use paramlist
      implicit none

      integer i,ina,ika
      real*8 dxyz(3),dd,rdmax,dtmp,dampmut
      real*8 sfamp,svamp,sprd,scosth,rfac
      real*8, allocatable:: vrrt(:,:)

      allocate(vrrt(3,natom))

      if(loopa.le.5) then
        dampmut=0.4d0
      else if(scosth.gt.0.7071d0) then
        dampmut=0.0d0
      else if(scosth.lt.0.0d0) then
        dampmut=1.0d0
      else
!c      dampmut=dampmu
        dampmut=0.4d0
      end if

      write(6,'(1x/1x,"damped quickmin  damping factor  mu=",f8.4)') 
     &  dampmut

! --- previous velocites
      do ina=1,natom
        if(iposfix(ina).eq.1) then
!          ika=katm(ina)
          rfac=tstep*omegai/(mass(ina)*mfac(ina))
          vrrt(:,ina)= rfac*matmul(frc(:,ina),bv)
          if(loopa.ge.3) then
            vrr(:,ina)=vrr(:,ina)-vrrt(:,ina)
          end if
        end if
      end do

! --- velocity*force for the whole system
      sfamp=0.0d0
      svamp=0.0d0
      sprd=0.0d0
      do ina=1,natom
        if(iposfix(ina).eq.1) then
!          ika=katm(ina)
          sfamp=sfamp+sum(frc(:,ina)**2)
          svamp=svamp+sum(matmul(uv,vrr(:,ina))**2)
          sprd=sprd+sum(matmul(uv,vrr(:,ina))*frc(:,ina))
        end if
      end do
      sfamp=dsqrt(sfamp)
      svamp=dsqrt(svamp)

      if(sfamp.gt.1.0d-12) then
        dtmp=sprd/sfamp**2
      else
        dtmp=0.0d0
      end if

      write(6,'(1x,"atom    dx          dy          dz"/
     & 1x,"------------------------------------------")')

      do ina=1,natom
      if(iposfix(ina).ne.0) then
!        ika=katm(ina)

! --- make the velocities parallel to the present forces
        if(sfamp.gt.1.0d-12) then
          vrr(:,ina)= matmul(frc(:,ina),bv)*omegai*sprd/(sfamp)**2
        else
          vrr(:,ina)= 0.0d0
        end if

        vrr(:,ina)= (vrrt(:,ina)+vrr(:,ina)*(1.0d0-dampmut))
     & /(1.0d0+dampmut)
        dxyz(:)=tstep*matmul(uv,vrr(:,ina))
        dd=dsqrt(sum(dxyz**2))
        if(dd.gt.rdmax) then
          dxyz=dxyz*rdmax/dd
          dd=rdmax
        end if

        ra(:,ina)=ra(:,ina)+dxyz(:)
        write(6,'(2x,i3,2x,3d12.4)') ina,dxyz

      else  ! iposfix
         rah(:,ina,imod2)=rah(:,ina,imod)
         vrr(:,ina)=0.0d0
         gradh(:,ina,imod)=0.0d0
      end if  ! iposfix

      end do

      write(6,'(1x,"------------------------------------------"/)')

      deallocate(vrrt)

      return
      end

! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**
      subroutine damped_md(rdmax,scosth)
! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**
      use paramlist
      implicit none

      integer i,ina,ika,itmp
      real*8 dxyz(3),dd,rdmax,scosth,dampmut,rfac
      real*8, allocatable:: vrrt(:,:)

      allocate(vrrt(3,natom))

      if(loopa.le.5) then
        dampmut=0.4d0
      else if(scosth.gt.0.7071d0) then
        dampmut=0.0d0
      else if(scosth.lt.0.0d0) then
        dampmut=1.0d0
      else
        dampmut=0.4d0
      end if

      write(6,'(1x/1x,"damped_md  damping factor  mu=",f8.4)') dampmut

      write(6,'(1x,"atom    dx          dy          dz"/
     & 1x,"------------------------------------------")')

      do ina=1,natom

      if(iposfix(ina).eq.1) then

!        ika=katm(ina)
! --- velocity at previous step
        rfac=tstep*omegai/(mass(ina)*mfac(ina))
        vrrt(:,ina)=rfac*matmul(frc(:,ina),bv)
        if(loopa.ge.3) then
          vrr(:,ina)=vrr(:,ina)-vrrt(:,ina)
        end if

        vrr(:,ina)= (vrrt(:,ina)+vrr(:,ina)*(1.0d0-dampmut))
     & /(1.0d0+dampmut)

        dxyz=tstep*matmul(uv,vrr(:,ina))
        dd=dsqrt(sum(dxyz**2))
        if(dd.gt.rdmax) then
          dxyz=dxyz*rdmax/dd
        end if

        ra(:,ina)=ra(:,ina)+dxyz(:)
        write(6,'(2x,i3,2x,4d12.4)') ina,dxyz
      else  ! iposfix
        rah(:,ina,imod2)=rah(:,ina,imod)
        vrr(i,ina)=0.0d0
        gradh(:,ina,imod)=0.0d0
      end if  ! iposfix

      end do

      write(6,'(1x,"------------------------------------------"/)')

      deallocate(vrrt)

      return
      end

! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**
!  2008.01.09  T.Tamura                                                *
!                                                                      *
      subroutine gdiis(mlin,rdmax,scosth,gnorm1m,gnorm0,gnorm,dtoter)
! *                                                                    *
! *    Geometry Optimization by DIIS                                   *
! *        "P.Csaszar and P.Pulay, J.Mole.Struct.114(1984)31-34"       *
! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**
      use paramlist
      implicit none

      integer i,ina,ika,mlin,mlint,ierror,iflag_md
      real*8 rdmax,aw
      real*8 scosth,prd,famp,vamp,gnorm,gnorm0,prdt,dtmp,dtoter,famp0

!c    integer imod1m
      real*8 gnorm1m

!      iuphess=0
      iuphess=1

!c    mlin=5

      write(6,'(1x,"imod=",i3,3x,"igdiis=",i5)') imod,igdiis

! 2008.02.27
      call gdistance(mlin,mlint)

! --- start of GDIIS

! 2008.03.03
        iflag_md=0
        if(igdiis.gt.1000) then
          if(mlint.ge.4.and.gnorm.lt.gnorm0) then
            go to 130
          else
            igdiis=0
            go to 110
          end if
        else if(scosth.lt.0.0d0) then

          if(mlint.ge.4) then
            go to 130
          else
            igdiis=0
            call resetvel
            go to 110
          end if
        else
          go to 110
        end if

! ----- simple-md

  110 call simplemd2(rdmax,iflag_md)

      igdiis=igdiis+1

      go to 1000

! ----- Relaxation by GDIIS

  130 ierror=0
      call gdiis0(mlin,mlint,ierror,rdmax)

      if (ierror.ne.0) then
          write(6,'(1x,"fault in gdiis")')
! 2008.02.28
        igdiis=1
        call resetvel
        go to 110
      end if

! 2008.02.27
      if(igdiis.lt.999) then
        igdiis=1001
      else
        igdiis=igdiis+1
      end if

      do ina=1,natom
        ra(:,ina)=rah(:,ina,imod2)
      end do

 1000 continue

      call outtstep

      return
      end

! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**
!  2008.01.09
      subroutine gdiis0(mlin,mlint,ierror,rdmax)
! *                                                                    *
! *     Main part of geometyr DIIS method                              *
! *       P.Csaszar and P.Pulay, J.Mole.Struct.114(1984)31-34          *
! *       O.Farkas and H.B.Schlegel, Phys.Chem.Chem.Phys.,4(2002)11-15 *
! *     pos1: Interpolated geometry                                    *
! *     grad1: Interpolated gradient                                   *
! *     err: Error matrix {e(i)=-H^(-1)*g(i)}                          *
! *     bb: Scalar product of the error vectors {BB(ij)=<e(i)|e(j)>}   *
! *     acc: Coefficients of Linear Conbination {c(i)}                 *
! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**
      use paramlist
      implicit none

      integer i,j,k,ina,mlin,info,ierror,mlint,mlin0
      integer ix,iy,iz,jj,il,ika
      integer, allocatable::indx(:)
      real*8 sumc,sabc,dtmp,aw,
     &  dxyz(3),dd,rwork1(1),rdmax,famp
      real*8, allocatable:: pos1(:,:),grad1(:,:),
     &  bb(:,:),acc(:),ipiv(:),rwork2(:),acc2(:,:)
      real*8, allocatable:: hgrad1(:,:)

      real*8 dtmp1,dtmp2,max_sabc,max_sabc_2
      integer mlintt
      integer gdiis_history

      gdiis_history=3

      max_sabc=20.0d0
!c    max_sabc=40.0d0
!c    max_sabc=100.0d0

      max_sabc_2=100.0d0

      write(6,'(1x/1x,"Geometry optimization by GDIIS  mlint=",i3)') 
     &  mlint

      if (iuphess.eq.1) then
        call bfgs(mlin)
      end if

      allocate(pos1(3,natom),grad1(3,natom),acc2(mlin,mlin))
      allocate(hgrad1(3,natom))

      acc2=0.0d0

      ierror=0
!c    mlin0=mlin
      mlin0=mlint

      write(6,'(1x,"coefficients of linear combination")')

  100 if(mlint.le.1) then
        go to 190
      end if

! ***** calculation of coefficients (acc,acc2)

      allocate(indx(mlint),bb(mlint,mlint),acc(mlint),
     &  ipiv(mlint),rwork2(mlint))

      do i=1,mlint
        j=mod((imod+mlin-mlint+i-1),mlin)+1
        indx(i)=j
      end do

      if (mlint.gt.mlin) then
        deallocate(indx,bb,acc,ipiv,rwork2)
        mlint=mlint-1
        mlin0=mlint
        go to 100
      end if

      do j=1,mlint
      do i=1,j
        bb(i,j)=0.0d0
      end do
      end do

      do ina=1,natom
        if(iposfix(ina).ge.1) then
          do j=1,mlint
          do i=1,j
        if(iuphess.eq.0) then
          do k=1,3
            bb(i,j)=bb(i,j)+gradh(k,ina,indx(i))*gradh(k,ina,indx(j))
          end do
        else if(iuphess.eq.1) then
          do k=1,3
            jj=k+(ina-1)*3
            bb(i,j)=bb(i,j)
     &  +gradh(k,ina,indx(i))*gradh(k,ina,indx(j))*(hessi2(jj,jj)**2)
          end do
        end if
          end do
          end do
        end if
      end do

      dtmp=1.0d0/bb(1,1)
      do j=1,mlint
      do i=1,j
        bb(i,j)=bb(i,j)*dtmp
      end do
      end do

      call dsytrf('U',mlint,bb,mlint,ipiv,rwork1,1,info)
      if (info.ne.0) then
        write(6,*) 'error in dsytrf'
        ierror=ierror+1
        deallocate(indx,bb,acc,ipiv,rwork2)
        mlint=mlint-1
        go to 1000
      end if

      call dsytri('U',mlint,bb,mlint,ipiv,rwork2,info)
      if (info.ne.0) then
        write(6,*) 'error in dsytri'
        ierror=ierror+1
        deallocate(indx,bb,acc,ipiv,rwork2)
        mlint=mlint-1
        go to 1000
      end if

      do j=1,mlint
      do i=1,j
        bb(j,i)=bb(i,j)
      end do
      end do

      do i=1,mlint
        acc(i)=0.d0
        do j=1,mlint
        acc(i)=acc(i)+bb(i,j)
        end do
      end do

      sumc=0.d0
      do j=1,mlint
      do i=1,mlint
        sumc=sumc+bb(i,j)
      end do
      end do

      do i=1,mlint
        acc(i)=acc(i)/sumc
      end do

      sabc=0.0d0
      do i=1,mlint
        sabc=sabc+dabs(acc(i))
      end do
      write(6,'(3x,"dimension=",i3,3x,"sabc=",d12.4)') mlint,sabc

      do i=1,mlint
        acc2(i,mlint)=acc(i)
      end do

      deallocate(indx,bb,acc,ipiv,rwork2)
      mlint=mlint-1
      go to 100

  190 mlintt=mlin0
  191 if(mlintt.le.(gdiis_history-1)) then
        write(6,'(1x,"fault in sabc")')
        ierror=1
        go to 1000
      end if
      sabc=0.0d0
      do i=1,mlintt
        sabc=sabc+dabs(acc2(i,mlintt))
      end do
      if((mlintt.ge.4.and.sabc.gt.max_sabc).or.
     &  (mlintt.le.3.and.sabc.gt.max_sabc_2)) then
        mlintt=mlintt-1
        go to 191
      end if

! ********** relaxation

  192 write(6,621) 
  621 format(1x,'atom    dx          dy          dz          dd'/
     & 1x,'--------------------------------------------------------')

      do 200 ina=1,natom

      if (iposfix(ina).ge.1) then  ! iposfix

      mlint=mlin0

  201 if(mlint.le.1) then
        go to 230
      end if

      if(mlint.gt.mlin0) then
        mlint=mlint-1
        go to 201
      end if

      sabc=0.0d0
      do i=1,mlint
        sabc=sabc+dabs(acc2(i,mlint))
      end do
      if((mlint.ge.4.and.sabc.gt.max_sabc).or.
     &  (mlint.le.3.and.sabc.gt.max_sabc_2)) then
        mlint=mlint-1
        go to 201
      end if

! --- Interpolated geometry and gradient

      allocate(indx(mlint))
      do i=1,mlint
        j=mod((imod+mlin-mlint+i-1),mlin)+1
        indx(i)=j
      end do
      do k=1,3
        pos1(k,ina)=0.0d0
        grad1(k,ina)=0.0d0
      end do
      do k=1,3
      do i=1,mlint
        pos1(k,ina)=pos1(k,ina)+acc2(i,mlint)*rah(k,ina,indx(i))
        grad1(k,ina)=grad1(k,ina)+acc2(i,mlint)*gradh(k,ina,indx(i))
      end do
      end do
      deallocate(indx)

      if(iuphess.eq.0) then
        hgrad1(1,ina)=grad1(1,ina)
        hgrad1(2,ina)=grad1(2,ina)
        hgrad1(3,ina)=grad1(3,ina)
      else if(iuphess.eq.1) then
        ix=1+(ina-1)*3
        iy=2+(ina-1)*3
        iz=3+(ina-1)*3
        hgrad1(1,ina)=hessi2(ix,ix)*grad1(1,ina)
        hgrad1(2,ina)=hessi2(iy,iy)*grad1(2,ina)
        hgrad1(3,ina)=hessi2(iz,iz)*grad1(3,ina)
      end if
        hgrad1(1,ina)=2.0d0*grad1(1,ina)
        hgrad1(2,ina)=2.0d0*grad1(2,ina)
        hgrad1(3,ina)=2.0d0*grad1(3,ina)

      dxyz(:)=pos1(:,ina)-hgrad1(:,ina)-rah(:,ina,imod)

      dd=dsqrt(sum(dxyz**2))
      if(dd.gt.rdmax) then
        dxyz=dxyz*rdmax/dd
        dd=rdmax
      end if
      vrr(:,ina)=0.0d0
      go to 210

! --- simple relaxation

  230 dxyz=-gradh(:,ina,imod)

      dd=dsqrt(sum(dxyz(:)**2))
      if(dd.gt.rdmax) then
        dxyz=dxyz*rdmax/dd
        dd=rdmax
      end if
      vrr(:,ina)=0.0d0

  210 rah(:,ina,imod2)=rah(:,ina,imod)+dxyz(:)

      write(6,'(2x,i3,2x,4d12.4,2x,i2)') 
     &  ina,dxyz(1:3),dd,mlint

      else   ! iposfix

      rah(:,ina,imod2)=rah(:,ina,imod)
      vrr(:,ina)=0.0d0

      end if  ! iposfix

  200 continue

      write(6,623)
  623 format(1x,'------------------------------',
     & '--------------------------')

      ierror=0

 1000 continue

      deallocate(pos1,grad1,acc2)
      deallocate(hgrad1)

      return
      end

! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**
!c    subroutine simplemd2(rdmax,scosth,dtoter)
      subroutine simplemd2(rdmax,iflag_md)
!       with acceleration
! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**
      use paramlist
      implicit none

      integer ina,i,ika,iflag_md

      real*8 dxyz(3),dd,rdmax,dtmp,rfac

      if(igdiis.ge.11) then
        dtmp=1.2d0**(igdiis/5-1)
      else
        dtmp=1.0d0
      end if


      write(6,'(1x/1x,"simple MD  damping factor =",f8.4)') dtmp

      write(6,'(1x,"atom    dx          dy          dz"/
     & 1x,"------------------------------------------")')

      do ina=1,natom
!        ika=katm(ina)

        if(iposfix(ina).ne.0) then

          rfac=tstep*omegai/(mass(ina)*mfac(ina))
          vrr(:,ina)=vrr(:,ina)-rfac*matmul(frc(:,ina),bv)
          vrr(:,ina)=vrr(:,ina)*dtmp
          vrr(:,ina)=vrr(:,ina)+rfac*matmul(frc(:,ina),bv)

          dxyz=tstep*matmul(uv,vrr(:,ina))
          dd=dsqrt(sum(dxyz**2))
          if(dd.gt.rdmax) then
            dxyz=dxyz*rdmax/dd
            dd=rdmax
            vrr(:,ina)=vrr(:,ina)*rdmax/dd
          end if
           ra(:,ina)=ra(:,ina)+dxyz

          write(6,'(2x,i3,2x,3d12.4)') ina,dxyz

        else
          rah(:,ina,imod2)=rah(:,ina,imod)
          vrr(:,ina)=0.0d0
          gradh(:,ina,imod)=0.0d0
        end if

      end do

      write(6,'(1x,"------------------------------------------"/)')

      return

      end

! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**
      subroutine bfgs(mlin)
! *                                                                    *
! *     Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm              *
! *       "Numerical Recipes in C (in japanese) p315"                  *
! *        update using the last min0(mlin,4) geometries               *
! *         see F.Eckert et al., J.Comp.Chem.18,1473(1997)             *
! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**
      use paramlist
      implicit none

      integer ina,il,i,j,mlin,mlint,ii,jj,ndim,ix,iy,iz
      integer, allocatable:: indx(:),indx0(:)
      real*8 fac,fae,fad,hesmax,hesmin
      real*8, allocatable:: ddp(:),ddg(:),hdg(:)

      hesmin=0.1d0
      hesmax=10.0d0

      ndim=3*natom

      allocate(ddp(ndim),ddg(ndim),hdg(ndim))

      write(6,'(1x,"update of the inverse Hessian")')

      do ii=1,ndim
      do jj=1,ndim
        hessi2(jj,ii)=0.0d0
      end do
        hessi2(ii,ii)=1.0d0
      end do

      if(loopa.le.2) then
        go to 1000
      end if
      mlint=5

      allocate(indx(mlint),indx0(mlint))

      do i=1,mlint
        j=mod((imod+mlin-mlint+i-1),mlin)+1
        indx(i)=j
        indx0(i)=mod((j+mlin-2),mlin)+1
      end do

      do 100 i=1,mlint

      if (i.eq.1) go to 100

        if((loopa.le.mlin).and.(i.eq.2)) then
          go to 100
        end if

        do ina=1,natom
        do il=1,3
          ii=il+(ina-1)*3
          ddp(ii)=rah(il,ina,indx(i))-rah(il,ina,indx0(i))
          ddg(ii)=gradh(il,ina,indx(i))-gradh(il,ina,indx0(i))
        end do
        end do

        do ii=1,ndim
          hdg(ii)=0.0d0
        do jj=1,ndim
          hdg(ii)=hdg(ii)+hessi2(ii,jj)*ddg(jj)
        end do
        end do

        fac=0.0d0
        fae=0.0d0
        do ii=1,ndim
          fac=fac+ddg(ii)*ddp(ii)
          fae=fae+ddg(ii)*hdg(ii)
        end do

        fac=1.0d0/fac
        fad=1.0d0/fae

        do ii=1,ndim
          ddg(ii)=fac*ddp(ii)-fad*hdg(ii)
        end do

        do ii=1,ndim
        do jj=1,ndim
          hessi2(ii,jj)=hessi2(ii,jj)+fac*ddp(ii)*ddp(jj)
     &    -fad*hdg(ii)*hdg(jj)+fae*ddg(ii)*ddg(jj)
        end do
        end do

  100 continue

      deallocate(indx,indx0)

 1000 continue

      do ii=1,ndim
        hessi2(ii,ii)=dmin1(hessi2(ii,ii),hesmax)
        hessi2(ii,ii)=dmax1(hessi2(ii,ii),hesmin)
      end do

      do ina=1,natom
        ix=1+(ina-1)*3
        iy=2+(ina-1)*3
        iz=3+(ina-1)*3
      end do

      deallocate(ddp,ddg,hdg)

      return
      end

! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**
!  2008.02.01  T.Tamura
      subroutine resetvel
! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**

      use paramlist
      implicit none

      integer ina,ika

      write(6,'(1x,"velocity is resetted to zero")')

      do ina=1,natom
!        ika=katm(ina)
        vrr(:,ina)=tstep*matmul(frc(:,ina),bv)
     & *omegai/(mass(ina)*mfac(ina))
      end do

      return
      end


! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**
!  2008.02.27  T.Tamura
      subroutine gdistance(mlin,mlint)
!       eliminate any geometris more than a certain distance
!      from the current geometry. (dafalut = 0.2 a.u.)
! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**

      use paramlist
      implicit none

      integer i,ina,imodt,mlin,mlint
      real*8 ddmax,dxyz(3),dd

      write(6,'(1x,"distance from current geometry")')

      mlint=0

      do i=1,(loopa-2)
        imodt=mod((loopa-1-i),mlin)+1
        ddmax=0.0d0
        do ina=1,natom
          if(iposfix(ina).ne.0) then
          dxyz=rah(:,ina,imodt)-rah(:,ina,imod)
          dd=dsqrt(sum(dxyz**2))
          ddmax=dmax1(ddmax,dd)
          end if
        end do
        if(ddmax.lt.0.2d0) then
          mlint=i
        end if
      end do
      write(6,'(3x,"mlint=",i4)') mlint

      return
      end
! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**
      subroutine outtstep
      use paramlist
      implicit none

      open(7,file='tstep.txt',status='unknown')
      write(7,*) 'tstep'
      write(7,*) tstep
      write(7,*) 'igdiis'
      write(7,*) igdiis
      write(7,*) 'lattice_mass'
      write(7,*) mcell
      write(7,*) 'md_gdiis'
      write(7,*) imdgdiis 
      write(7,'(1x,"damping_factor")') 
      write(7,*) dampmu
      close(7)

      return
      end
! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**
      subroutine intstep
      use paramlist
      implicit none

      open(7,file='tstep.txt',status='old')
      read(7,*)
      read(7,*) tstep
      read(7,*)
      read(7,*) igdiis
      read(7,*)
      read(7,*) mcell
      read(7,*)
      read(7,*) imdgdiis 
      read(7,*) 
      read(7,*) dampmu
      close(7)

      return
      end
! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**
