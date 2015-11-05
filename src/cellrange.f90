
module m_cellrange
  implicit none

  private:: normalize,f_cross_x
contains

  !２次元でいうと、四隅の点に対して、(p,q+0.5),(p+0.5,q)の平面への距離を計算する。
  function make_supercell(plat,rcut) result (supercell)
    implicit none
    real(8),intent(in):: plat(3,3)
    real(8),intent(in):: rcut
    integer:: supercell(3)  ! return value

    integer:: i1,i2,i3
    integer:: imax
    real(8):: p(3),corner1(3),corner2(3), corner3(3),nvec(3),r

    integer:: fid=6

    supercell=[0,0,0] ! initialize
    imax=10

    ! (p,q,r+0.5)などへの平面を定義

    write(fid,*)'(100)'
    ! (1,0,0) plane
    corner1=plat(:,2)
    corner2=plat(:,3)
    corner3=f_cross_x(corner1,corner2)
    nvec=normalize(corner3)
    i1=0;i2=0;i3=0
    do i1=1,imax 
       ! 平面上の点を定義
        p(:)=(i1*0.5)*plat(:,1)+i2*plat(:,2)+(i3)*plat(:,3)
!         (l n - p ). n =0
        r= abs(sum(nvec*p))
        write(fid,'(2(a,3F12.5,1x),2(a,F12.5,1x))') 'p=',p,'n=',nvec,'r=',r,'rcut=',rcut
	if (r>rcut) then
           supercell(1)=i1
	   exit
        endif
    enddo

    write(fid,*)'(010)'
    ! (0,1,0) plane
    corner1=plat(:,3)
    corner2=plat(:,1)
    corner3=f_cross_x(corner1,corner2)
    nvec=normalize(corner3)
    i1=0;i2=0;i3=0
    do i2=1,imax 
       ! 平面上の点を定義
        p(:)=i1*plat(:,1)+(i2*0.5d0)*plat(:,2)+(i3)*plat(:,3)
!         (l n - p ). n =0
        r= sum(nvec*p) 
        write(fid,'(2(a,3F12.5,1x),2(a,F12.5,1x))') 'p=',p,'n=',nvec,'r=',r,'rcut=',rcut
        if (r>rcut) then
           supercell(2)=i2
	   exit
        endif
    enddo 

    write(fid,*)'(001)'
    !(0,0,1)平面
    corner1=plat(:,1)
    corner2=plat(:,2)
    corner3=f_cross_x(corner1,corner2)
    nvec=normalize(corner3)
    i1=0;i2=0;i3=0
    do i3=1,imax
       ! 平面上の点を定義
        p(:)=i1*plat(:,1)+i2*plat(:,2)+(i3*0.5d0)*plat(:,3)
!	  (l n - p ). n =0
	r= sum(nvec*p)
        write(fid,'(2(a,3F12.5,1x),2(a,F12.5,1x))') 'p=',p,'n=',nvec,'r=',r,'rcut=',rcut
        if (r>rcut) then
           supercell(3)=i3
	   exit
        endif
    enddo

    supercell=supercell

  end function make_supercell


  function f_cross_x(a,b) result(c)
    implicit none
    real(8) :: a(3),b(3),c(3)
    c(1)=a(2)*b(3)-a(3)*b(2)
    c(2)=a(3)*b(1)-a(1)*b(3)
    c(3)=a(1)*b(2)-a(2)*b(1)
  end function f_cross_x

  function normalize(a) result(b)
    implicit none
    real(8) :: a(3),b(3)
    real(8):: r
    r=sum(a**2)
    r=sqrt(r)
    b=a/r
  end function normalize

end module m_cellrange

#ifdef MAIN
!---------------------
program range
  use m_cellrange
  implicit none

  real(8):: plat(3,3),rcut,superlat(3,3)
  integer:: icell(3),i

  rcut=2.1d0
  plat(:,1)= [ 1.0d0,1.0d0,0.0d0]
  plat(:,2)= [0.0d0, 1.0d0, 0.0d0]
  plat(:,3)= [0.0d0,0.0d0,3.0d0] 
  write(*,*)'rcut=',rcut
  write(*,*)'plat='
  write(*,'(3F15.5)')plat
  icell=  make_supercell(plat,rcut)
  write(*,*)'isupercell=',icell
  do i=1,3
  superlat(:,i)= icell(i)*plat(:,i)
  enddo
  write(*,*)'supercell='
  write(*,'(3F15.5)')superlat
end program range

#endif
