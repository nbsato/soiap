module rfc5
  use paramlist
  use cif_module
  implicit none

  ! parameters for optimization
  logical, parameter :: fix_first_atom = .false. ! fix coords of the first atom for debug
  logical, parameter :: fix_atoms = .false.  ! fix atoms coords for debug
  logical, parameter :: fix_units = .false.  ! fix unit vectors for debug
  integer, parameter :: gdiis_start_step = 2 ! start of GDIIS steps
  integer, parameter :: gdiis_steps = 6      ! depth of GDIIS steps
  real(8), parameter :: criterion_max_diff = 0.10d0 ! modulation criterion for diff[Bohr]

  ! module arrays for optimization history
  real(8), allocatable :: coord_recent(:,:,:)
  real(8), allocatable :: force_recent(:,:,:)
  real(8), allocatable :: Hessian(:,:,:,:)

contains
  ! optimize atom coordinates and cell vectors
  subroutine updt_coordcell_RFC5( loop )
    integer, intent(in) :: loop ! optimization iteration

    integer :: size ! total size of positions
    integer :: i, j, k, ia
    integer :: info
    real(8) :: force_mean(3), max_diff

    ! local variables for GDIIS
    real(8), allocatable :: overlap_ff(:,:), beta(:)
    integer, allocatable :: ipiv(:)
    real(8), allocatable :: work(:), eigen(:)
    real(8), allocatable :: force_gdiis(:,:)
    real(8), allocatable :: coord_gdiis(:,:)

    ! local variables for BFGS
    real(8), allocatable :: coord_diff(:,:)
    real(8), allocatable :: force_diff(:,:)
    real(8), allocatable :: direc_diff(:,:)
    real(8), allocatable :: coord_bfgs(:,:)
    real(8) :: prod_rf, prod_rd
    real(8), allocatable :: RF(:,:,:,:)
    real(8), allocatable :: Z(:,:,:)
    integer              :: M
    integer, allocatable :: iwork(:), ifail(:)

    !--- OptSW variables
    ! number of atom: QMD%natom
    ! atom coord: QMD%ra(:,ia) (abosolute)
    ! atom coord: QMD%rr(:,ia) (relative)
    ! atom force: QMD%frc(:,ia)
    !
    ! cell vectors: QMD%uv(:,1), QMD%uv(:,2), QMD%uv(:,3)
    ! recp vectors: QMD%bv(:,1), QMD%bv(:,2), QMD%bv(:,3)
    ! cell volume: QMD%omega
    ! cell stress: QMD%strs(:,:)
    !---

    size = QMD%natom + 3 ! sum of atoms and vectors

    ! allocate arrays for recording history
    if( loop == 1 ) then
       allocate( coord_recent(3,size,gdiis_steps) )
       allocate( force_recent(3,size,gdiis_steps) )
       coord_recent = 0.0d0
       force_recent = 0.0d0

       allocate( Hessian(3,size,3,size) )
       Hessian = 0.0d0
    end if

    ! allocate local arrays
    allocate( force_gdiis(3,size), coord_gdiis(3,size) )
    allocate( force_diff(3,size), coord_diff(3,size) )
    allocate( direc_diff(3,size), coord_bfgs(3,size) )

    !---- GDIIS stage

    ! shift the history of coords and forces
    do i=gdiis_steps-1, 1, -1
       coord_recent(:,:,i+1) = coord_recent(:,:,i)
       force_recent(:,:,i+1) = force_recent(:,:,i)
    end do

    ! record the latest coords and forces of atoms
    do ia=1, QMD%natom
       coord_recent(:,ia,1) = QMD%ra(:,ia) ! use absolute coordinate
       force_recent(:,ia,1) = QMD%frc(:,ia) ! use force in absolute coordinate
    end do

    ! remove numerical errors in the sum of force_recent(:,:,1)
    force_mean(:) = 0.0d0
    do ia=1, QMD%natom
       force_mean(:) = force_mean(:) + force_recent(:,ia,1)
    end do
    force_mean(:) = force_mean(:)*(1.0d0/QMD%natom)
    do ia=1, QMD%natom
       force_recent(:,ia,1) = force_recent(:,ia,1) - force_mean(:)
    end do

    ! record the latest coords and forces of unit vectors
    do ia=1, 3
       coord_recent(:,QMD%natom+ia,1) = QMD%uv(:,ia) ! use absolute coordinate
       force_recent(:,QMD%natom+ia,1) = QMD%omega * matmul( QMD%strs(:,:), QMD%bv(:,ia) )
    end do

    if( fix_atoms ) then
       do ia=1, QMD%natom
          force_recent(:,ia,1) = 0.0d0
       end do
    else if( fix_first_atom ) then
       ia=1
       force_recent(:,ia,1) = 0.0d0
    end if

    if( fix_units ) then
       do ia=1, 3
          force_recent(:,QMD%natom+ia,1) = 0.0d0
       end do
    end if

    if( loop<gdiis_start_step ) then
       coord_gdiis(:,:) = coord_recent(:,:,1)
       force_gdiis(:,:) = force_recent(:,:,1)
    else
       allocate( overlap_ff(gdiis_steps,gdiis_steps), beta(gdiis_steps) )
       ! calc matrix of l.h.s of the linear equation
       do j=1, gdiis_steps
          do i=j, gdiis_steps
             if( i>loop .or. j>loop ) then
                if( i== j ) then
                   overlap_ff(i,j) = 1.0d0
                else
                   overlap_ff(i,j) = 0.0d0
                end if
             else
                overlap_ff(i,j) = sum( force_recent(:,:,i)*force_recent(:,:,j) )
             end if
          end do
       end do

       ! calc vector of r.h.s of the linear equation
       do i=1, gdiis_steps
          if( i>loop ) then
             beta(i) = 0.0d0
          else
             beta(i) = 1.0d0
          end if
       end do

       ! solve the linear equation
       allocate( ipiv(gdiis_steps), work(gdiis_steps*4) )
       call dsysv( 'L', gdiis_steps, 1, overlap_ff, gdiis_steps, ipiv, &
            beta, gdiis_steps, work, 4*gdiis_steps, info )
       deallocate( ipiv, work )
       beta(:) = beta(:)*(1.d0/sum(beta(:)))

       ! calc GDIIS trial coords and forces as a linear combination
       coord_gdiis(:,:) = 0.0d0
       force_gdiis(:,:) = 0.0d0
       do i=1, gdiis_steps
          coord_gdiis(:,:) = coord_gdiis(:,:) + beta(i)*coord_recent(:,:,i)
          force_gdiis(:,:) = force_gdiis(:,:) + beta(i)*force_recent(:,:,i)
       end do

       deallocate( overlap_ff, beta )
    end if

    !---- BFGS stage
    if( loop == 1 ) then
       Hessian = 0.0d0
       do ia=1, size
          do k=1, 3
             Hessian(k,ia,k,ia) = 1.0d0
          end do
       end do
    else
       ! calc difference of coords and forces from the previous step
       coord_diff(:,:) = coord_recent(:,:,1) - coord_recent(:,:,2)
       force_diff(:,:) = force_recent(:,:,1) - force_recent(:,:,2)

       ! operate the hessian on the difference of the forces
       do ia=1, size
          do k=1, 3
             direc_diff(k,ia) = sum(Hessian(k,ia,:,:)*coord_diff(:,:))
          end do
       end do

       ! calc inner products
       prod_rf = sum(coord_diff(:,:)*force_diff(:,:))
       prod_rd = sum(coord_diff(:,:)*direc_diff(:,:))

       ! update the hessian
       do ia=1, size
          do k=1, 3
             Hessian(k,ia,:,:) = Hessian(k,ia,:,:) &
                  - (force_diff(k,ia)/prod_rf) * force_diff(:,:) &
                  - (direc_diff(k,ia)/prod_rd) * direc_diff(:,:)
          end do
       end do
    end if

    !---- RF stage

    ! calc matrix of RF
    ! [memo] RF matrix is a little larger than 3size+1 x 3size+1.
    ! extra elements in this matrix are not accessed in dsyev.
    allocate( RF(3,size+1,3,size+1), Z(3,size+1,1) )
    RF(:,1:size,:,1:size) = Hessian(:,1:size,:,1:size)
    RF(1,size+1,:,1:size) = (-1.0d0)*force_gdiis(:,1:size)
    RF(1,size+1,1,size+1) = 0.0d0

    ! solve the lowest eigen vector
    allocate( work(8*(3*size+1)), eigen(3*size+1) )
    allocate( iwork(5*(3*size+1)), ifail(3*size+1) )

    call dsyevx( 'V', 'I', 'L', 3*size+1, RF, 3*(size+1), &
         0.0d0, 0.0d0, 1, 1, 0.0d0, M, &
         eigen, Z, 3*(size+1), work, 8*(3*size+1), iwork, ifail, info )

    ! normalize eigen vector as ( vec{coord_bfgs}, 1 )
    coord_bfgs(:,:) = Z(:,1:size,1)*(1.0d0/Z(1,size+1,1))

    deallocate( work, eigen )
    deallocate( iwork, ifail )
    deallocate( RF, Z )

    ! calc diff of coord of unit vectors on the next step
    do ia=1, 3
       coord_diff(:,QMD%natom+ia) = &
            coord_gdiis(:,QMD%natom+ia) + coord_bfgs(:,QMD%natom+ia) &
            - QMD%uv(:,ia)
    end do

    ! calc diff of coord of atoms on the next step
    do ia=1, QMD%natom
       coord_diff(:,ia) = coord_gdiis(:,ia) + coord_bfgs(:,ia) - QMD%ra(:,ia)
    end do

    ! symmetrize diff of coord of atoms
    if( CIF_canSymmetrize() ) then
       ! translate from absolute to relative coordinates
       do ia=1, QMD%natom
          coord_diff(:,ia) = matmul(coord_diff(:,ia),QMD%bv(:,:))
       end do
       call CIF_symmetrizeDirection( coord_diff )
       ! translate from relative to absolute coordinates
       do ia=1, QMD%natom
          coord_diff(:,ia) = matmul(QMD%uv(:,:),coord_diff(:,ia))
       end do
    end if

    ! find maximum coord_diff of unit vectors and atoms in absolute coordinates
    max_diff = maxval(abs(coord_diff(:,:)))

    ! modulation coord_diff if it is too large
    if( max_diff > criterion_max_diff ) then
       coord_diff(:,:) = coord_diff(:,:)*(criterion_max_diff/max_diff)
    end if

    ! update coords of unit vectors on the next step
    do ia=1, 3
       if( .not. fix_units ) then
          QMD%uv(:,ia) = QMD%uv(:,ia) + coord_diff(:,QMD%natom+ia)
       end if
    end do

    ! update reciprocal vectors
    call cross_x(QMD%uv(:,1),QMD%uv(:,2),QMD%bv(:,3))
    call cross_x(QMD%uv(:,2),QMD%uv(:,3),QMD%bv(:,1))
    call cross_x(QMD%uv(:,3),QMD%uv(:,1),QMD%bv(:,2))
    QMD%omega=dot_product(QMD%bv(:,3),QMD%uv(:,3))
    QMD%omegai=1.d0/QMD%omega
    QMD%bv=QMD%bv*QMD%omegai

    ! update coords of atoms on the next step
    do ia=1, QMD%natom
       if( fix_atoms .or. (fix_first_atom.and.ia==1) ) then
          if( .not. fix_units ) then
             QMD%ra(:,ia) = matmul(QMD%uv(:,:),QMD%rr(:,ia))
          end if
       else
          QMD%ra(:,ia) = matmul(QMD%uv(:,:),QMD%rr(:,ia)) + coord_diff(:,ia)
          QMD%rr(:,ia) = matmul(QMD%ra(:,ia),QMD%bv(:,:))
       end if
    end do

    ! clear local arrays
    deallocate( force_diff, coord_diff )
    deallocate( direc_diff, coord_bfgs )
    deallocate( force_gdiis, coord_gdiis )

  end subroutine updt_coordcell_RFC5
end module rfc5
