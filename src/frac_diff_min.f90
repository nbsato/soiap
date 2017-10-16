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

subroutine frac_diff_min(p,p1,diff)
   implicit none
   real(8),intent(in):: p(3)
   real(8),intent(in)::p1(3)
   integer:: i,i1,i2,i3,ilist(3)
   real(8):: diffmin(3),diff(3),difflen,difflendir
   integer::direction(3),ilistdir(3)

   diff=p-p1
   direction=0
   do i=1,3
   if (diff(i)>0.5d0) direction(i)=1
   if (diff(i)<-0.5d0) direction(i)=-1
   enddo

#define NEWDIFF 1

#if !defined(NEWDIFF)
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
#endif

   diff=p-p1
   difflen=sum(diff**2)
   ilistdir=[0,0,0]
      diff= p-p1-direction
      if (difflen>sum(diff**2)) then
        difflen=sum(diff**2)
        ilistdir=direction
      endif
     
#if defined(NEWDIFF)
    ilist=ilistdir
#else
   if (ilistdir(1)/=ilist(1) .or. &
   ilistdir(2)/=ilist(2) .or. &
   ilistdir(3)/=ilist(3) ) then
      diff=p-p1-ilist; difflen=sum(diff**2)
      diff=p-p1-ilistdir; difflendir=sum(diff**2)

   if (abs(difflen-difflendir)> 1.0d-6) then 
       write(6,*)'failed'
      write(6,*)'p=',p
      write(6,*)'p1=',p1
      write(6,*)'direction=',direction,p-p1
      write(6,*)  'ilistdir/=ilist'
      diff=p-p1-ilist; difflen=sum(diff**2)
      write(6,*) 'ilist=',ilist,difflen
      diff=p-p1-ilistdir; difflendir=sum(diff**2)
      write(6,*) 'ilistdir=',ilistdir,difflendir

   do i3=-1,1
   do i2=-1,1
   do i1=-1,1
     diff = p-p1-[ i1,i2,i3]
     difflen= sum(diff**2)
     write(6,'(3i3,4F10.5)')i1,i2,i3,difflen,diff
   enddo
   enddo
   enddo

   stop "abort"

   endif

   endif
#endif

   i1=ilist(1); i2=ilist(2); i3=ilist(3)
#if 1
   diff = p - p1 - [ real(i1,kind=8),real(i2,kind=8),real(i3,kind=8) ]
#else
   p1= p1+ [ real(i1,kind=8),real(i2,kind=8),real(i3,kind=8) ]
   diff= p-p1
#endif

end subroutine frac_diff_min

