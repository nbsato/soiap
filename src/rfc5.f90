module rfc5
  use paramlist
  use cif_module
  implicit none

  ! depth of GDIIS history
  integer, parameter   :: gdiis_steps = 6

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

    ! record the latest coords and forces of atoms and unit vectors
    do ia=1, QMD%natom
       coord_recent(:,ia,1) = QMD%rr(:,ia) ! use relative coordinate
       force_recent(:,ia,1) = QMD%frc(:,ia)
    end do
    do ia=1, 3
       coord_recent(:,QMD%natom+ia,1) = QMD%uv(:,ia)
       force_recent(:,QMD%natom+ia,1) = QMD%omega * matmul( QMD%strs(:,:), QMD%bv(:,ia) )
    end do

    if( loop<gdiis_steps ) then
       coord_gdiis(:,:) = coord_recent(:,:,1)
       force_gdiis(:,:) = force_recent(:,:,1)
    else
       allocate( overlap_ff(gdiis_steps,gdiis_steps), beta(gdiis_steps) )
       ! calc matrix of l.h.s of the linear equation
       do j=1, gdiis_steps
          do i=j, gdiis_steps
             overlap_ff(i,j) = sum( force_recent(:,:,i)*force_recent(:,:,j) )
          end do
       end do

       ! calc vector of r.h.s of the linear equation
       do i=1, gdiis_steps
          beta(i) = 1.0d0
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
    allocate( work(8*(3*size+1)), eigen(3*(size+1)) )
    allocate( iwork(5*(3*size+1)), ifail(3*size+1) )
    call dsyevx( 'V', 'I', 'L', 3*size+1, RF, 3*(size+1), &
         0.0d0, 0.0d0, 1, 1, 0.0d0, M, &
         eigen, Z, 3*(size+1), work, 8*(3*size+1), iwork, ifail, info )

    ! normalize eigen vector as ( vec{coord_bfgs}, 1 )
    coord_bfgs(:,:) = Z(:,1:size,1)*(1.0d0/Z(1,size+1,1))

    deallocate( work, eigen )
    deallocate( iwork, ifail )
    deallocate( RF, Z )

    ! update coords of unit vectros on the next step
    do ia=1, 3
       QMD%uv(:,ia) = coord_gdiis(:,QMD%natom+ia) + coord_bfgs(:,QMD%natom+ia)
    end do

    ! calc diff of coord of atoms on the next step
    do ia=1, QMD%natom
       coord_diff(:,ia) = coord_gdiis(:,ia) + coord_bfgs(:,ia) - QMD%rr(:,ia)
    end do

    ! symmetrize diff of coord of atoms
    if( CIF_canSymmetrize() ) then
       call CIF_symmetrizeDirection( coord_diff )
    end if

    ! update coords of atoms on the next step
    do ia=1, QMD%natom
       QMD%rr(:,ia) = QMD%rr(:,ia) + coord_diff(:,ia)
    end do

    ! update related variables
    do ia=1, QMD%natom
       QMD%ra(:,ia) = matmul(QMD%uv(:,:),QMD%rr(:,ia))
    end do
    call cross_x(QMD%uv(:,1),QMD%uv(:,2),QMD%bv(:,3))
    call cross_x(QMD%uv(:,2),QMD%uv(:,3),QMD%bv(:,1))
    call cross_x(QMD%uv(:,3),QMD%uv(:,1),QMD%bv(:,2))
    QMD%omega=dot_product(QMD%bv(:,3),QMD%uv(:,3))
    QMD%omegai=1.d0/QMD%omega
    QMD%bv=QMD%bv*QMD%omegai

    ! clear local arrays
    deallocate( force_diff, coord_diff )
    deallocate( direc_diff, coord_bfgs )
    deallocate( force_gdiis, coord_gdiis )

  end subroutine updt_coordcell_RFC5
end module rfc5
