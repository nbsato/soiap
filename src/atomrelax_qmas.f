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
        stop 'imd=1: under debug'
!        call simple_md
      else if(iabs(imd).eq.2) then
        stop 'imd=2: under debug'
!        if(loopa.le.3) then
!          call quickmin(rdmax,scosth)
!        else
!          call gdiis(mlin,rdmax,scosth,gnorm1m,gnorm0,gnorm,dtoter)
!        end if
      else if(iabs(imd).eq.3) then
        call simple_relax
      else if(iabs(imd).eq.4) then
        call quench_md(rdmax)
      else if(iabs(imd).eq.5) then
        stop 'imd=5: under debug'
!        if(loopa.le.3) then
!          call quickmin(rdmax,scosth)
!        else
!          call damped_md(rdmax,scosth)
!        end if
      else if(iabs(imd).eq.6) then
        stop 'imd=6: under debug'
!        call quickmin(rdmax,scosth)
!      else
!        call simple_relax
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
