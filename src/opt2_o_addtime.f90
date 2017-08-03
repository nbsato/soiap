program main
  use m_showtime
  use paramlist
  use frcfield
  use rfc5
  implicit none

  integer i,j,ia,ika
  integer:: loopa,loopc
  real*8:: rfac
  type(t_showtime)::  timer

  open(901,file='log.struc',form='formatted')
  open(902,file='log.tote',form='formatted')
  write(902,*)'# loopc,loopa,tote,fmax,tote/natom,omega,omega/natom'
  open(903,file='log.frc',form='formatted')
  open(904,file='log.strs',form='formatted')

  call input

  if (QMD%imd/=0.and.abs(QMD%imd)/=3.and.abs(QMD%imd)/=4) then
     write(*,*)'md_mode should be 0, 3 or 4'
     write(*,*)'md_mode=0: FIRE mode'
     write(*,*)'md_mode=3: simple relax'
     write(*,*)'md_mode=4: quenched MD (default)'
     stop 'QMD%imd error: other modes in structure_opt are under debug.'
  endif
  if (QMD%imdc/=0.and.QMD%imdc/=2.and.QMD%imdc/=3) then
     write(*,*)'md_mode_cell should be 0 2 or 3'
     write(*,*)'md_mode_cell=0: FIRE'
     write(*,*)'md_mode_cell=2: quenched MD (default)'
     write(*,*)'md_mode_cell=3: RFC5'
     stop 'QMD%imdc error: other modes in structure_opt are under development.'
  endif

  if( QMD%imdc == 3 ) then
     do loopc=1,QMD%nloopc
        QMD%loopc=loopc
        do loopa=1,QMD%nloopa
           QMD%loopa=loopa
           write(6,*)'QMD%loopc,QMD%loopa=',QMD%loopc,QMD%loopa,' out of ', &
                QMD%nloopc,QMD%nloopa

           call timer%start()
           call tote_frc_strs
           call timer%stop()
           call timer%show("tote_frc_strs:")

           do i=1,3
              QMD%strs(i,i)=QMD%strs(i,i)-QMD%extstrs(i)
           enddo

           call output_tote(902)

           if (QMD%fmax<QMD%fth) then
              write(*,'("QMD%frc converged. QMD%loopc, QMD%loopa, QMD%fmax",2i5,e12.4)') &
                   QMD%loopc,QMD%loopa,QMD%fmax
              exit
           endif

           if (loopa/=QMD%nloopa) then
              call updt_coordcell_RFC5((QMD%loopc-1)*QMD%nloopa+QMD%loopa,.true.,.false.)
           endif
        enddo

        call strs_max

        call output_struc(901)
        call output_frc(903)
        call output_strs(904)

        if (QMD%smax<QMD%sth) then
           write(*,'("QMD%strs converged. QMD%loopc, QMD%smax",i5,e12.4)') &
                 QMD%loopc,QMD%smax
           if (QMD%fmax<QMD%fth) exit
        endif

        call updt_coordcell_RFC5((QMD%loopc-1)*QMD%nloopa+QMD%loopa,.true.,.true.)
     enddo
  else
     call tote_frc_strs
     QMD%uvo(:,:,1)=QMD%uv
     QMD%strso(:,:,1)=QMD%strs
     QMD%vuv=0.d0
     do loopc=1,QMD%nloopc ! relax unit cell
        QMD%loopc=loopc

        do loopa=1,QMD%nloopa ! relax internal coordinates
           QMD%loopa=loopa
           write(6,*)'QMD%loopc,QMD%loopa=',QMD%loopc,QMD%loopa,' out of ', &
                QMD%nloopc,QMD%nloopa

           call timer%start()
           call tote_frc_strs
           call timer%stop()
           call timer%show("tote_frc_strs:")

           do i=1,3
              QMD%strs(i,i)=QMD%strs(i,i)-QMD%extstrs(i)
           enddo

           ! debug: >
           call output_tote(902)
           if (QMD%fmax<QMD%fth) then 
              write(*,'("QMD%frc converged. QMD%loopc, QMD%loopa, QMD%fmax",2i5,e12.4)') &
                   QMD%loopc,QMD%loopa,QMD%fmax
              exit
           endif

           call updt_velocity_atom
           call updt_coord
        enddo ! QMD%loopa

        call strs_max

        call output_struc(901)
        call output_frc(903)
        call output_strs(904)
        if (QMD%smax<QMD%sth) then
           write(*,'("QMD%strs converged. QMD%loopc, QMD%smax",i5,e12.4)')&
                QMD%loopc,QMD%smax
           if (QMD%fmax<QMD%fth) exit
        endif

        call updt_velocity_cell
        call updt_cell
     enddo ! QMD%loopc
  end if ! imdc /= 3

  close(901)
  close(902)
  close(903)
  close(904)

  write(*,*)'*** QMD%loopc =',QMD%loopc
  write(*,'(a,10F20.10)')'tote,fmax,tote/natom,omega,omega/natom=',&
       QMD%tote,QMD%fmax,&
       QMD%tote/QMD%natom,QMD%omega*((bohr)**3),QMD%omega*(bohr**3)/QMD%natom 

  write(6,*)'!!!! normally end !!!!'

contains

  subroutine tote_frc_strs
    real*8 faemp

    if (QMD%ifrcf==1) then ! Stillinger-Weber for Si
       call Stillinger_Weber
    elseif (QMD%ifrcf==2) then ! Tsuneyuki potential for Si-O
       stop 'Tsuneyuki potential is under debug'
       call Tsuneyuki
    elseif (QMD%ifrcf==3) then ! ZRL potential for Si-O
       call ZRL
    elseif (QMD%ifrcf==4) then ! ADP for Nd-Fe-B
       call ADP_KWU14
    else   
       stop 'etot_frc_strs error'
    endif

    QMD%fmax=0.d0
    do i=1,QMD%natom
       faemp=dsqrt(sum(QMD%frc(:,i)**2))
       if (QMD%iposfix(i)/=0.and.faemp>QMD%fmax) then
          QMD%fmax=faemp
          QMD%iamax=i
       endif
       !
       !     if (QMD%iposfix(i)/=0.and.faemp>QMD%fcut) then
       !        QMD%frc(:,i)=QMD%frc(:,i)*QMD%fcut/faemp
       !     endif   
    enddo

    !  write(*,'(a,10F20.10)')'tote,fmax,tote/natom,omega,omega/natom=',&
    !QMD%tote,QMD%fmax,&
    !QMD%tote/QMD%natom,QMD%omega*((bohr)**3),QMD%omega*(bohr**3)/QMD%natom 

  end subroutine tote_frc_strs

  subroutine strs_max
    integer :: i,j

    QMD%smax=0.d0
    do i=1,3
       do j=1,3
          if (i==j) then
             QMD%smax=QMD%smax+QMD%strs(i,j)**2
          else   
             QMD%smax=QMD%smax+QMD%strs(i,j)**2*2.d0
          endif
       enddo
    enddo

    QMD%smax=dsqrt(QMD%smax)

  end subroutine strs_max

  subroutine updt_cell

    integer :: ifin,ret,ilength,icolumn,iatompos

    QMD%uvo(:,:,2)=QMD%uvo(:,:,1)
    QMD%uvo(:,:,1)=QMD%uv
    QMD%strso(:,:,2)=QMD%strso(:,:,1)
    QMD%strso(:,:,1)=QMD%strs

    !  QMD%vuv=QMD%vuv+0.5d0*QMD%tstep*matmul(QMD%strs,QMD%uv)/QMD%mcell

    if (QMD%imdc==0) then ! fire
       call lattice_fire(0.1d0)
    elseif (QMD%imdc==1) then ! steepest descent
       call latticerelax_sd
    elseif (QMD%imdc==2) then ! quenched MD
       call latticerelax
    else   
       call lattice_simple_relax
    endif

    call cross_x(QMD%uv(:,1),QMD%uv(:,2),QMD%bv(:,3))
    call cross_x(QMD%uv(:,2),QMD%uv(:,3),QMD%bv(:,1))
    call cross_x(QMD%uv(:,3),QMD%uv(:,1),QMD%bv(:,2))
    QMD%omega=dot_product(QMD%bv(:,3),QMD%uv(:,3))
    write(*,*)'QMD%omega =',QMD%omega
    QMD%omegai=1.d0/QMD%omega
    QMD%bv=QMD%bv*QMD%omegai

  end subroutine updt_cell

  subroutine updt_velocity_cell
    if (QMD%loopc>1) then
       !    QMD%vuv=QMD%vuv+QMD%tstep*matmul(QMD%strs,QMD%uv)/QMD%mcell
       QMD%vuv=QMD%vuv+QMD%tstepc*matmul(QMD%strs,QMD%uv)/QMD%mcell
    else
       !    QMD%vuv=QMD%vuv+QMD%tste*matmul(QMD%strs,QMD%uv)/QMD%mcell*0.5d0
       QMD%vuv=QMD%vuv+QMD%tstepc*matmul(QMD%strs,QMD%uv)/QMD%mcell*0.5d0
    endif
  end subroutine updt_velocity_cell

  subroutine updt_coord

    integer :: ia
    real*8 :: rfac,p1,p2,vrrmax,vrrabs
    real*8,allocatable :: ratmp(:,:,:)

    call atomrelax

  end subroutine updt_coord

  subroutine updt_velocity_atom
    integer :: ia
    real(8) :: rfac
    do ia=1,QMD%natom
       if (QMD%loopa>1) then
          rfac=QMD%tstep/(QMD%mass(ia)*QMD%mfac(ia))
       else
          rfac=QMD%tstep/(QMD%mass(ia)*QMD%mfac(ia))*0.5d0
       endif
       QMD%vrr(:,ia)=QMD%vrr(:,ia)+rfac*matmul(QMD%frc(:,ia),QMD%bv(:,:))
    enddo
  end subroutine updt_velocity_atom

  subroutine output_struc(ifo)
    integer :: ifo,i,j

    write(ifo,*)'*** unit vectors [Bohr]: QMD%loopc =',QMD%loopc
    do i=1,3
       write(ifo,"(3f23.16)")(QMD%uv(i,j),j=1,3)
    enddo
    write(ifo,*)'*** internal lattice coordinates: QMD%loopc =',QMD%loopc
    do j=1,QMD%natom
       write(ifo,"(3f23.16)")(QMD%rr(i,j),i=1,3)
    enddo

    return
  end subroutine output_struc

  subroutine output_tote(ifo)
    integer :: ifo,i,j

    write(ifo,'(2i5,5F20.10)')QMD%loopc,QMD%loopa,QMD%tote,QMD%fmax,&
         QMD%tote/QMD%natom,QMD%omega*((bohr)**3),QMD%omega*(bohr**3)/QMD%natom 

    !  write(ifo,'(a,2i5,5F20.10)')&
    !'loopa,loopc,tote,fmax,tote/natom,omega,omega/natom=',&
    !QMD%loopc,QMD%loopa,QMD%tote,QMD%fmax,&
    !QMD%tote/QMD%natom,QMD%omega*((bohr)**3),QMD%omega*(bohr**3)/QMD%natom 

    return
  end subroutine output_tote

  subroutine output_frc(ifo)
    integer :: ifo,i,j

    write(ifo,*)'*** atomic forces: QMD%loopc =',QMD%loopc
    do j=1,QMD%natom
       write(ifo,"(3f23.16)")(QMD%frc(i,j),i=1,3)
    enddo

    return
  end subroutine output_frc

  subroutine output_strs(ifo)
    integer :: ifo,i,j

    write(ifo,*)'QMD%strs'
    write(ifo,'("QMD%loopc, QMD%smax",i5,2e16.6,e23.10)')QMD%loopc,&
         QMD%smax,QMD%sth,QMD%tote
    do i=1,3
       write(ifo,"(3f20.10)")(QMD%strs(i,j),j=1,3)
    enddo

    return
  end subroutine output_strs

end program main
