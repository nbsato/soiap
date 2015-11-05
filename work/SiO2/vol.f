      implicit real*8(a-h,o-z)
      real*8 delta
      parameter (delta=0.01)
      a0=8.312678
      c0=5.332292

      do i=-5,5
         fac=1.0d0+delta*dble(i)
         a=a0*fac
         c=c0*fac
         vol=a*a*c/16
         write(*,*)'*** scale =',fac,vol
         write(*,900)a,0,0
         write(*,900)0,a,0
         write(*,900)0,0,c
         write(*,*)
      enddo   

 900  format("( ",3f16.12," )")

      end
