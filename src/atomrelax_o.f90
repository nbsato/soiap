      subroutine atomrelax
! QMD%imd  =0: FIRE
!          =3: simple ralaxation
!          =4: quenched MD                               
      use paramlist
      implicit none

      integer mlin,n,isp,il,ina
      real*8 rdmax,scosth,gnorm1m,gnorm0,gnorm,dtoter
      real*8, allocatable:: rrsave(:,:)

      allocate (rrsave(3,QMD%natom))
      rrsave=QMD%rr

! parameters
      mlin=QMD%nloopa
!      rdmax=0.25d0 ! read from inputfile. see subroutine input for default

! setup
      if ((QMD%loopa==1).and.(QMD%loopc==1)) then
        allocate (QMD%rah(3,QMD%natom,mlin),QMD%gradh(3,QMD%natom,mlin))
      end if
      if (QMD%loopa==1) QMD%npstv=0

! QMD%imod: current step,  QMD%imod0: previous step,  QMD%imod2: next step
      QMD%imod=mod((QMD%loopa-1),mlin)+1
      if (QMD%loopa>1) then
        QMD%imod0=mod((QMD%loopa-2),mlin)+1
      else
        QMD%imod0=mlin
      end if
      QMD%imod2=mod(QMD%loopa,mlin)+1

! relative to absolute
      QMD%ra=matmul(QMD%uv,QMD%rr)

! Save Positions and Graidents at the previous relaxation step
      do ina=1,QMD%natom
        QMD%rah(:,ina,QMD%imod)=QMD%ra(:,ina)
        QMD%gradh(:,ina,QMD%imod)=-QMD%frc(:,ina)
      end do

! atomic relaxation
      if(iabs(QMD%imd).eq.3) then
        call simple_relax
      else if(iabs(QMD%imd).eq.4) then
        call quench_md(QMD%rdmax)
      else
        call fire(QMD%rdmax)
      end if

! absolute to relative
      do ina=1,QMD%natom
         QMD%rr(:,ina)=matmul(QMD%ra(:,ina),QMD%bv)
      enddo

      QMD%toter0=QMD%tote

! for fixed atoms
      do ina=1,QMD%natom
        if (QMD%iposfix(ina).eq.0) then
          QMD%rr(:,ina)=rrsave(:,ina)
        end if
      end do
      QMD%ra=matmul(QMD%uv,QMD%rr)
      deallocate (rrsave)

      return
      end subroutine atomrelax
! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**
      subroutine simple_md
      use paramlist
      implicit none
      integer ina
      real*8 dra(3)

      do ina=1,QMD%natom
        dra=QMD%tstep*matmul(QMD%uv,QMD%vrr(:,ina))
        QMD%ra(:,ina)=QMD%ra(:,ina)+dra(:)
!        write(6,'(i3,2x,3d14.6)') ina,dra(1:3)
      end do

      return
      end subroutine simple_md
! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**
      subroutine simple_relax
      use paramlist
      implicit none
      integer ina
      real*8 dra(3)

      if((QMD%loopc.eq.1).and.(QMD%loopa.eq.1)) return

      do ina=1,QMD%natom
        if(QMD%iposfix(ina).eq.1) then
        dra=QMD%tstep*matmul(QMD%uv,QMD%vrr(:,ina))
        QMD%ra(:,ina)=QMD%ra(:,ina)+dra(:)
!        write(6,'(i3,2x,3d14.6)') ina,dra(1:3)
        end if
        QMD%vrr(:,ina)=0.0d0
      end do

      return
      end subroutine simple_relax
! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**
      subroutine quench_md(rdmax)
      use paramlist
      implicit none
      integer i,ina
      real*8 dra(3),dd,rdmax,prd

      if((QMD%loopc.eq.1).and.(QMD%loopa.eq.1)) return

      do ina=1,QMD%natom
        if(QMD%iposfix(ina).eq.1) then
          prd=sum(matmul(QMD%uv,QMD%vrr(:,ina))*QMD%frc(:,ina))
          if (prd.lt.(-tol)) then
             QMD%rah(:,ina,QMD%imod2)=QMD%rah(:,ina,QMD%imod)
             QMD%vrr(:,ina)=0.0d0
          end if

          dra(:)=QMD%tstep*matmul(QMD%uv,QMD%vrr(:,ina))
          dd=dsqrt(sum(dra**2))
          if(dd.gt.rdmax) then
            dra=dra*rdmax/dd
            dd=rdmax
          end if

          QMD%ra(:,ina)=QMD%ra(:,ina)+dra(:)
!          write(6,'(i3,2x,3d14.6)') ina,dxyz(1:3)
        end if
      end do

      return
      end subroutine quench_md
! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**
      subroutine fire(rdmax)
! Bitzek et al., PRL97, 170201 (2006)
! F1: calculate P = F.v
! F2: set v -> (1-alpha)v + alpha F^ |v|, where ^ denoting a unit vector
! F3: if P>0 and the number of steps since P was negative is larger thab Nmin,
!    increase the time step dt -> min(dt f_inc, dtmax) and decrease 
!    alpha -> alpha f_alpha
! F4: if P<=0, decrease time step dt -> dt f_dec, freeze the system v->0
!    and set alpha back to alpha_start
! F5: return to MD
!
! Parameters
! Nmin=5, f_inc=1.1, f_dec=0.5, alpha_start=0.1, f_alpha=0.99
! dt_max ~ 10 dt_MD
      use paramlist
      implicit none
      integer i,ina
      real*8 p,vnorm,fnorm,va(3,QMD%natom),dra(3),dd,rdmax
! parameters
      integer nmin
      real*8 finc,fdec,alp0,falp,dtmax
!
      if (QMD%loopc==1.and.QMD%loopa==1) then
        QMD%fire_nmin=5
        QMD%fire_finc=1.1d0
        QMD%fire_fdec=0.5d0
        QMD%fire_alp0=0.1d0
        QMD%fire_falp=0.99d0
        QMD%fire_dtmax=QMD%tstep0*10.0d0
        QMD%fire_alp=QMD%fire_alp0
        return
      endif
! F1
! p=F.v
      va=0.0d0
      p=0.0d0
      vnorm=0.0d0
      fnorm=0.0d0
      do ina=1,QMD%natom
        if(QMD%iposfix(ina).eq.1) then
           va(:,ina)=matmul(QMD%uv,QMD%vrr(:,ina))
           p=p+sum(va(:,ina)*QMD%frc(:,ina))
           vnorm=vnorm+sum(va(:,ina)**2)
           fnorm=fnorm+sum(QMD%frc(:,ina)**2)
        endif   
      enddo  
      if (p.gt.0.0d0) then 
         QMD%npstv=QMD%npstv+1
      else
         QMD%npstv=0
      endif   
      vnorm=sqrt(vnorm)
      fnorm=sqrt(fnorm)
! F2
      do ina=1,QMD%natom
        if(QMD%iposfix(ina).eq.1) then
           va(:,ina)=(1.0d0-QMD%fire_alp)*va(:,ina)+ &
           QMD%fire_alp*QMD%frc(:,ina)/fnorm*vnorm
        else   
           va(:,ina)=0.0d0
        endif
      enddo  
! F3
      if ((p.gt.0.0d0).and.(QMD%npstv.gt.QMD%fire_nmin)) then
         QMD%tstep=min(QMD%tstep*QMD%fire_finc,QMD%fire_dtmax)
         QMD%fire_alp=QMD%fire_alp*QMD%fire_falp
      endif   
! F4
      if (p.le.0.0d0) then
         QMD%tstep=QMD%tstep*QMD%fire_fdec
         va(:,:)=0.0d0
         QMD%fire_alp=QMD%fire_alp0
      endif   
! MD
      do ina=1,QMD%natom
        if(QMD%iposfix(ina).eq.1) then
          dra(:)=QMD%tstep*matmul(QMD%uv,QMD%vrr(:,ina))
          dd=dsqrt(sum(dra**2))
          if(dd.gt.rdmax) then
            dra=dra*rdmax/dd
            dd=rdmax
          endif 
          QMD%ra(:,ina)=QMD%ra(:,ina)+dra(:)
!          write(6,'(i3,2x,3d14.6)') ina,dxyz(1:3)
        end if
        QMD%vrr(:,ina)=matmul(va(:,ina),QMD%bv(:,:))
      enddo   

! debug:
      write(*,'(a,i5,3f16.8)')'fire',QMD%loopa,p,vnorm/QMD%natom,fnorm/QMD%natom

      return
      end subroutine fire
! **+****1****+****2****+****3****+****4****+****5****+****6****+****7**
