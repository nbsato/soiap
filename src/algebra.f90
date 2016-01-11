subroutine ddijdr(dij,rij,ddijdri,ddijdrj)
  implicit none
  real(8) :: dij,rij(3),ddijdri(3),ddijdrj(3)
  ddijdri=-rij/dij
  ddijdrj= rij/dij
end subroutine ddijdr

subroutine dcosjikdr(rij,rik,dij,dik,cosjik,dcdri,dcdrj,dcdrk)
  implicit none
  real(8) :: rij(3),rik(3),dij,dik,cosjik,dcdri(3),dcdrj(3),dcdrk(3)
  real(8) :: eij(3),eik(3)
  eij=rij/dij
  eik=rik/dik
  dcdrj=(eik-eij*cosjik)/dij
  dcdrk=(eij-eik*cosjik)/dik
  dcdri=-(dcdrj+dcdrk)
end subroutine dcosjikdr

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
