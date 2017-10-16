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

module m_showtime

type t_showtime
integer:: t1,t2,t_rate,t_max
contains
  procedure:: start => start_timer
  procedure:: stop=> end_timer
  procedure:: show=>show_time
end type

private start_timer, end_timer, show_time

contains

subroutine start_timer(v)
implicit none
class (t_showtime):: v
call system_clock(v%t1)
end subroutine start_timer

subroutine end_timer(v)
implicit none
class (t_showtime):: v
call system_clock(v%t2,v%t_rate,v%t_max)
end subroutine end_timer


subroutine show_time(v,s) !t1,t2,t_rate,t_max)
implicit none
class (t_showtime):: v
character(*):: s
!integer:: t1,t2,t_rate,t_max
integer::diff
  if ( v%t2 < v%t1 ) then
    diff = (v%t_max - v%t1) + v%t2 + 1
  else
    diff = v%t2 - v%t1
  endif
  write(6, "(2A,F10.3,A)") "show_time: ",trim(s), diff/dble(v%t_rate)," sec"
end subroutine show_time

real(8) function rtime(v) !t1,t2,t_rate,t_max)
implicit none
type(t_showtime):: v
!integer:: t1,t2,t_rate,t_max
integer::diff
  if ( v%t2 < v%t1 ) then
    diff = (v%t_max - v%t1) + v%t2 + 1
  else
    diff = v%t2 - v%t1
  endif
  rtime= diff/dble(v%t_rate) 
end function rtime

end module m_showtime
