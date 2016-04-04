!!!!! subroutines for relaxing atomic positions 

! for the time being, use atomrelax_qmas_o.f


!!!!! subroutines for optimizing the unit cell
subroutine latticerelax
  use paramlist

  real*8 :: auv(3,3,2),dtmp,c1,c2
  integer :: i,j

  do i=1,2
    auv(:,:,i)=matmul(QMD%strso(:,:,i),QMD%uvo(:,:,i))
  enddo

  do j=1,3
  do i=1,3
     if ((auv(i,j,1)*auv(i,j,2))<(-tol*1.d-10)) then
        if (QMD%loopc>1) then
!           dtmp=-auv(i,j,2)*QMD%uvo(i,j,1)/(auv(i,j,1)-auv(i,j,2)) &
!                +auv(i,j,1)*QMD%uvo(i,j,2)/(auv(i,j,1)-auv(i,j,2)) &
!                -QMD%uv(i,j)
           c1=-auv(i,j,2)/(auv(i,j,1)-auv(i,j,2))
           c1=dmin1(c1,1.d0)
           c1=dmax1(c1,0.d0)
           c2=1.d0-c1
           dtmp=c1*QMD%uvo(i,j,1)+c2*QMD%uvo(i,j,2)-QMD%uv(i,j)
        else   
!           dtmp=0.5d0*QMD%tstep*QMD%tstep*auv(i,j,1)/QMD%mcell
           dtmp=0.5d0*QMD%tstepc*QMD%tstepc*auv(i,j,1)/QMD%mcell
        endif
        dtmp=dmax1(dtmp,-0.1d0)
        dtmp=dmin1(dtmp, 0.1d0)
        QMD%uv(i,j)=QMD%uv(i,j)+dtmp
        QMD%vuv(i,j)=0.d0
     else
!        QMD%uv(i,j)=QMD%uv(i,j)+QMD%tstep*QMD%vuv(i,j)
        QMD%uv(i,j)=QMD%uv(i,j)+QMD%tstepc*QMD%vuv(i,j)
     endif   
  enddo ! i   
  enddo ! j   

end subroutine latticerelax

subroutine lattice_fire(rdmax)
  use paramlist
  implicit none
  integer i
  real*8 p,vnorm,fnorm,duv(3),dd,rdmax,f(3)
! parameters
      integer nmin
      real*8 finc,fdec,alp0,falp,dtmax
!
      if (QMD%loopc==1) then
        QMD%fire_nminc=5
        QMD%fire_fincc=1.1d0
        QMD%fire_fdecc=0.5d0
        QMD%fire_alp0c=0.1d0
        QMD%fire_falpc=0.99d0
        QMD%fire_dtmaxc=QMD%tstep0*10.0d0
        QMD%fire_alpc=QMD%fire_alp0c
        return
      endif
! F1
! p=F.v
! F=matmul(strs*uv)
      p=0.0d0
      vnorm=0.0d0
      fnorm=0.0d0
      do i=1,3
        f=matmul(QMD%strs,QMD%uv(:,i))
        p=p+sum(QMD%vuv(:,i)*f)
        vnorm=vnorm+sum(QMD%vuv(:,i)**2)
        fnorm=fnorm+sum(f**2)
      enddo  
      if (p.gt.0.0d0) then 
         QMD%npstvc=QMD%npstvc+1
      else
         QMD%npstvc=0
      endif   
      vnorm=sqrt(vnorm)
      fnorm=sqrt(fnorm)
! F2
      do i=1,3
        QMD%vuv(:,i)=(1.0d0-QMD%fire_alpc)*QMD%vuv(:,i)+ &
        QMD%fire_alpc*f(:)/fnorm*vnorm
      enddo  
! F3
      if ((p.gt.0.0d0).and.(QMD%npstv.gt.QMD%fire_nminc)) then
         QMD%tstepc=min(QMD%tstepc*QMD%fire_fincc,QMD%fire_dtmaxc)
         QMD%fire_alpc=QMD%fire_alpc*QMD%fire_falpc
      endif   
! F4
      if (p.le.0.0d0) then
         QMD%tstepc=QMD%tstepc*QMD%fire_fdecc
         QMD%vuv=0.0d0
         QMD%fire_alpc=QMD%fire_alp0c
      endif   
! MD
      do i=1,3
          duv(:)=QMD%tstepc*QMD%vuv(:,i)
          dd=dsqrt(sum(duv**2))
          if(dd.gt.rdmax) then
            duv=duv*rdmax/dd
            dd=rdmax
          endif 
          QMD%uv(:,i)=QMD%uv(:,i)+duv(:)
!          write(6,'(i3,2x,3d14.6)') ina,duv(1:3)
      enddo   

! debug:
      write(*,'(a,i5,3f16.8)')'lattice_fire',QMD%loopc,p,vnorm/3,fnorm/3

end subroutine lattice_fire

subroutine lattice_simple_relax
  use paramlist
! QMD%uv=QMD%uv+QMD%tstep*QMD%vuv
 QMD%uv=QMD%uv+QMD%tstepc*QMD%vuv
 QMD%vuv=0.d0
end subroutine lattice_simple_relax

subroutine latticerelax_sd
! modified from latticerelax_cg in QMAS
  use paramlist
  integer i,j,k
  real*8 lambt1,lambt2,lambu1,lambu2,lambda,dtmp1,dtmp2,gamma
  real*8 acoe,bcoe,dhmmax,lamb0,dd

!  trial step
  lambt1=25.0d0
  lambt2=50.0d0

!  maximum of lambda
  lambu1=100.0d0
  lambu2=10.0d0

  if(QMD%loopc==2) then
     QMD%isd_cell=0
     QMD%guv(:,:,0)=0.0d0
  endif

  QMD%guv(:,:,2)=matmul(QMD%strs,QMD%uv)

  if(QMD%loopc>=3) then
     dtmp1=0.0d0
     dtmp2=0.0d0
     do j=1,3
     do i=1,3
        dtmp1=dtmp1+QMD%guv(i,j,0)*QMD%guv(i,j,0)
        dtmp2=dtmp2+QMD%guv(i,j,0)*QMD%guv(i,j,2)
     end do
     end do
     write(6,'(1x,"g0*s_cell=",d16.8,2x,"g2*s_cell=",d16.8)') dtmp1,dtmp2
     if(dabs(dtmp2).lt.dabs(dtmp1)*0.01d0) QMD%isd_cell=0
  endif

  if(QMD%isd_cell==4) QMD%isd_cell=0

  QMD%isd_cell=QMD%isd_cell+1

! --- new search direction (sdv)

  if(QMD%isd_cell==1) then

     if(QMD%loopc<=3) then
        lambda=lambt1
     else
        lambda=lambt2
     endif

     dd=0.0d0
     do j=1,3
     do i=1,3
        dd=dmax1(dd,dabs(QMD%guv(i,j,2)))
     enddo
     enddo
     lambda=dmin1(lambda,(0.5d0/(dd*bohr)))

     QMD%uv(:,:)=QMD%uv(:,:)+lambda*QMD%guv(:,:,2)
     QMD%duv0(:,:)=lambda*QMD%guv(:,:,2)

     QMD%guv(:,:,0)=QMD%guv(:,:,2)
     QMD%guv(:,:,1)=QMD%guv(:,:,2)

! --- search the minimum along the search direction (sdv)

  else if(QMD%isd_cell >= 2) then

     dtmp1=0.0d0
     dtmp2=0.0d0
     do j=1,3
     do i=1,3
        dtmp1=dtmp1+QMD%guv(i,j,0)*QMD%guv(i,j,1)
        dtmp2=dtmp2+QMD%guv(i,j,0)*QMD%guv(i,j,2)
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

     QMD%uv(:,:)=QMD%uv(:,:)+lambda*QMD%duv0(:,:)
     QMD%duv0(:,:)=lambda*QMD%duv0(:,:)

     QMD%guv(:,:,1)=QMD%guv(:,:,2)

  end if ! QMD%isd_cell

  write(6,'(1x,"QMD%isd_cell=",i2,2x,"lambda=",d16.8)') QMD%isd_cell,lambda

  write(6,*)

  return
end subroutine latticerelax_sd
