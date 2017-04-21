!!----------------------------------------------------------------
!! functions for read CIF in OpenMX Viewer
!!----------------------------------------------------------------
module cif_module
  implicit none

  real(8), parameter :: M_PI       = 3.141592653589793238462d0 !< Pi, circular constant

  !!----------------
  !! structure for cell properties.
  !!----------------
  type Cell_type
     real(8) :: length_a
     real(8) :: length_b
     real(8) :: length_c
     real(8) :: angle_alpha
     real(8) :: angle_beta
     real(8) :: angle_gamma
  end type Cell_type

  !!----------------
  !! structure for symmetric operation
  !!----------------
  type Operation_type
     real(8) :: transM(3,3)
     real(8) :: transV(3)
  end type Operation_type

  !!----------------
  !! structure for atomic properties.
  !!----------------
  type Atom_type
     character(len=4) :: symbol
     real(8) :: coords(3)
  end type Atom_type

  !!----------------
  !! structure for CIF.
  !!----------------
  type CIF_type
     !! type of cell
     type(Cell_type) :: cell

     !! all symmetric operations as matrix and vector.
     type(Operation_type), pointer :: voperation(:) => null()

     !! all atoms recorded in CIF.
     type(Atom_type), pointer :: vatom_cif(:) => null()

     !! conventional unit vectors in cartesian coordinates.
     real(8) :: unit_conv(3,3)

     !! all atoms in the conventional unit cell.
     type(Atom_type), pointer :: vatom_conv(:) => null()

     !! symmetrize matrix for atom motion
     type(Operation_type), allocatable :: symmetrize(:,:)
  end type CIF_type
  type(CIF_type), save :: cif

contains
  !!----------------
  !! functions for element symbols.
  !!----------------
  function getAtomicNumber( symbol ) result(number)
    character(len=*), intent(in) :: symbol
    integer :: number

    select case( symbol )
    case( "H"); number= 1; case("He"); number= 2; case("Li"); number= 3;
    case("Be"); number= 4; case( "B"); number= 5; case( "C"); number= 6;
    case( "N"); number= 7; case( "O"); number= 8; case( "F"); number= 9;
    case("Ne"); number=10; case("Na"); number=11; case("Mg"); number=12;
    case("Al"); number=13; case("Si"); number=14; case( "P"); number=15;
    case( "S"); number=16; case("Cl"); number=17; case("Ar"); number=18;
    case( "K"); number=19; case("Ca"); number=20; case("Sc"); number=21;
    case("Ti"); number=22; case( "V"); number=23; case("Cr"); number=24;
    case("Mn"); number=25; case("Fe"); number=26; case("Co"); number=27;
    case("Ni"); number=28; case("Cu"); number=29; case("Zn"); number=30;
    case("Ga"); number=31; case("Ge"); number=32; case("As"); number=33;
    case("Se"); number=34; case("Br"); number=35; case("Kr"); number=36;
    case("Rb"); number=37; case("Sr"); number=38; case( "Y"); number=39;
    case("Zr"); number=40; case("Nb"); number=41; case("Mo"); number=42;
    case("Tc"); number=43; case("Ru"); number=44; case("Rh"); number=45;
    case("Pd"); number=46; case("Ag"); number=47; case("Cd"); number=48;
    case("In"); number=49; case("Sn"); number=50; case("Sb"); number=51;
    case("Te"); number=52; case( "I"); number=53; case("Xe"); number=54;
    case("Cs"); number=55; case("Ba"); number=56; case("La"); number=57;
    case("Ce"); number=58; case("Pr"); number=59; case("Nd"); number=60;
    case("Pm"); number=61; case("Sm"); number=62; case("Eu"); number=63;
    case("Gd"); number=64; case("Tb"); number=65; case("Dy"); number=66;
    case("Ho"); number=67; case("Er"); number=68; case("Tm"); number=69;
    case("Yb"); number=70; case("Lu"); number=71; case("Hf"); number=72;
    case("Ta"); number=73; case( "W"); number=74; case("Re"); number=75;
    case("Os"); number=76; case("Ir"); number=77; case("Pt"); number=78;
    case("Au"); number=79; case("Hg"); number=80; case("Tl"); number=81;
    case("Pb"); number=82; case("Bi"); number=83; case("Po"); number=84;
    case("At"); number=85; case("Rn"); number=86; case("Fr"); number=87;
    case("Ra"); number=88; case("Ac"); number=89; case("Th"); number=90;
    case("Pa"); number=91; case( "U"); number=92; case("Np"); number=93;
    case("Pu"); number=94; case("Am"); number=95; case("Cm"); number=96;
    case("Bk"); number=97; case("Cf"); number=98; case("Es"); number=99;
    case("Fm"); number=100; case("Md"); number=101; case("No"); number=102;
    case("Lr"); number=103;
    case default
       write(*,'(a,a)') '# Error!: unknown element symbol:', trim(symbol) 
       stop
    end select
  end function getAtomicNumber

  !!----------------
  !! functions for 3D coordinates.
  !!----------------
  function Position( x, y, z )
    real(8), intent(in) :: x, y, z
    real(8) :: Position(3)
  
    Position(:) = [ x, y, z ]
  end function Position

  function Position_outer_product( p1, p2 ) result(product)
    real(8), intent(in) :: p1(3), p2(3)
    real(8) :: product(3)

    product(:) = [ &
         p1(2)*p2(3) - p2(2)*p1(3), &
         p1(3)*p2(1) - p2(3)*p1(1), &
         p1(1)*p2(2) - p2(1)*p1(2) ]
  end function Position_outer_product

  function Position_inner_product( p1, p2 ) result(product)
    real(8), intent(in) :: p1(3), p2(3)
    real(8) :: product

    product = sum(p1(:)*p2(:))
  end function Position_inner_product

  function Position_match( p1, p2, eps )
    real(8), intent(in) :: p1(3)
    real(8), intent(in) :: p2(3)
    real(8), intent(in) :: eps

    logical :: Position_match
    real(8) :: dx, dy, dz
    
    dx = abs(p1(1)-p2(1))
    if( dx>eps .and. abs(dx-1.0d0)>eps ) then
       Position_match = .false.
       return
    end if

    dy = abs(p1(2)-p2(2))
    if( dy>eps .and. abs(dy-1.0d0)>eps ) then
       Position_match = .false.
       return
    end if

    dz = abs(p1(3)-p2(3))
    if( dz>eps .and. abs(dz-1.0d0)>eps ) then
       Position_match = .false.
       return
    end if

    Position_match = .true.
  end function Position_match


  !!----------------
  !! functions for atomic properties.
  !!----------------
  function Atom_match( atom1, atom2 )
    type(Atom_type), intent(in) :: atom1, atom2
    real(8), parameter :: eps = 2.0d-4
    logical :: Atom_match

    Atom_match = &
         atom1%symbol == atom2%symbol .and. &
         Position_match( atom1%coords, atom2%coords, eps )

  end function Atom_match

  function Atom_find( vatom, atom ) result(found)
    type(Atom_type), pointer    :: vatom(:)
    type(Atom_type), intent(in) ::  atom
    logical :: found
    
    integer :: n

    if( .not. associated(vatom) ) then
       found = .false.
       return
    end if

    do n=1, size(vatom)
       if( Atom_match(vatom(n),atom) ) then
          found = .true.
          return
       end if
    end do
    found = .false.
  end function Atom_find


  subroutine push_back_Atom( v, e )
    type(Atom_type), pointer :: v(:)
    type(Atom_type), intent(in) :: e

    integer :: N
    type(Atom_type), allocatable :: vnew(:)

    if( associated(v) ) then
       N = size(v)
    else
       N = 0
    end if
    allocate( vnew(N+1) )
    vnew(1:N) = v(1:N)
    vnew(N+1) = e

    if( associated(v) ) deallocate(v)
    allocate( v(N+1) )
    v(1:N+1) = vnew(1:N+1)
    deallocate( vnew )
  end subroutine push_back_Atom


  !!----------------
  !! functions for symmetric operation.
  !!----------------
  function Operation_construct( xyz ) result(operation)
    character(len=*), intent(in) :: xyz
    type(Operation_type) :: operation

    character(len=256) :: axis(3)
    character(len=256) :: fraction
    integer :: is, ie, n
    real(8) :: numerator, denominator

    is = 1
    if( xyz(is:is) == "'" ) then
       is = is + 1
    end if

    ie = index(xyz(is:),",")
    if( ie<=0 ) goto 1010
    axis(1) = xyz(is:is-1+ie-1)

    is = is-1+ie+1
    ie = index(xyz(is:),",")
    if( ie<=0 ) goto 1010
    axis(2) = xyz(is:is-1+ie-1)
    
    is = is-1+ie+1
    ie = index(xyz(is:),"'")
    if( ie<=0 ) then
       axis(3) = xyz(is:)
    else
       axis(3) = xyz(is:is-1+ie-1)
    end if

    do n=1, 3
       if( index(trim(axis(n)), "-x")>0 ) then
          operation%transM(n,1) = -1.0d0
       else if( index(trim(axis(n)), "x")>0 ) then
          operation%transM(n,1) = +1.0d0
       else
          operation%transM(n,1) =  0.0d0
       end if

       if( index(trim(axis(n)), "-y")>0 ) then
          operation%transM(n,2) = -1.0d0
       else if( index(trim(axis(n)), "y")>0 ) then
          operation%transM(n,2) = +1.0d0
       else
          operation%transM(n,2) =  0.0d0
       end if

       if( index(trim(axis(n)), "-z")>0 ) then
          operation%transM(n,3) = -1.0d0
       else if( index(trim(axis(n)), "z")>0 ) then
          operation%transM(n,3) = +1.0d0
       else
          operation%transM(n,3) =  0.0d0
       end if
    end do

    do n=1, 3
       is = verify(trim(axis(n)), "+-/.0123456789", back=.true.)
       fraction = axis(n)(is+1:)
       
       if( fraction == "" ) then
          operation%transV(n) = 0.0d0;
          cycle
       end if
       
       ie = verify(trim(fraction), "+-.0123456789")
       if( ie<=0 ) then
          read(fraction(1:),*) numerator
       else
          read(fraction(1:ie-1),*) numerator
       end if
       
       is = verify(trim(fraction), ".0123456789", back=.true.)
       if( is<=0 ) then
          denominator = 1.0d0
       else
          read(fraction(is+1:),*) denominator
       end if
       
       operation%transV(n) = numerator/denominator
       if( operation%transV(n) < 0.0d0 ) then
          operation%transV(n) = operation%transV(n) + 1.0d0
       end if
    end do
    return

1010 continue 
    write(*,*) "broken operation string", trim(xyz)
    stop
  end function Operation_construct

  function Operation_operate( operation, p_org ) result(p_new)
    type(Operation_type), intent(in) :: operation
    real(8), intent(in) :: p_org(3)
    real(8) :: p_new(3)

    integer :: n

    p_new(:) = matmul( operation%transM(:,:), p_org(:) )
    p_new(:) = p_new(:) + operation%transV(:)

    do n=1, 3
       do while( p_new(n)< 0.0d0 )
          p_new(n) = p_new(n) + 1.0d0
       end do
       do while( p_new(n)>=1.0d0 )
          p_new(n) = p_new(n) - 1.0d0
       end do
    end do
  end function Operation_operate


  subroutine push_back_Operation( v, e )
    type(Operation_type), pointer :: v(:)
    type(Operation_type), intent(in) :: e

    integer :: N
    type(Operation_type), allocatable :: vnew(:)

    if( associated(v) ) then
       N = size(v)
    else
       N = 0
    end if
    allocate( vnew(N+1) )
    vnew(1:N) = v(1:N)
    vnew(N+1) = e

    if( associated(v) ) deallocate(v)
    allocate( v(N+1) )
    v(1:N+1) = vnew(1:N+1)
    deallocate( vnew )
  end subroutine push_back_Operation

  !!----------------
  !! functions for processing line data.
  !!----------------
  function uncomment( line ) result(flag)
    character(len=*), intent(inout) :: line
    logical :: flag

    integer :: i
    integer :: len

    flag = .false.
    len = len_trim(line)

    if( len == 0 ) then
       flag = .true.
       return
    end if

    do i=1, len
       if( line(i:i) == '!' .or. line(i:i) == '#' .or. line(i:i) == '$' ) then
          line(i:len) = ' '
          if( i==1 ) then
             flag = .true.
          end if
          exit
       end if
       if( line(i:i) == '=' ) then
          line(i:i) = ' '
       end if
       if( line(i:i) == "'" .and. line(i+1:i+1) /= " " ) then
          line(i:i) = ' '
       end if
    end do

    if( len_trim(line) == 0 ) then
       flag = .true.
    end if

    return
  end function uncomment

  subroutine split( value, buf )
    character(len=*), intent(out) :: value(:)
    character(len=*), intent(in) :: buf

    integer :: lbuf, i, j, lval

    lbuf = len(buf)
    lval = 0
    j=1
    value(:) = ""
    do i=1, lbuf
       if( buf(i:i) == " " ) then
          if( lval>0 ) then
             j = j + 1
             lval = 0
          end if
       else
          value(j) = trim(value(j)) // buf(i:i)
          lval = lval + 1
       end if
    end do

  end subroutine split

  !!----------------
  !! functions for CIF.
  !!----------------

  !!----
  !! load CIF parameters.
  !!----
  subroutine CIF_readData( filename )
    character(len=*), intent(in) :: filename

    integer, parameter :: iunit = 101
    logical :: ex
    character(256) :: buf, tag, value(10)
    integer :: i

    logical :: loop_atom, loop_symm
    integer :: column, column_max
    integer :: column_xyz, column_x, column_y, column_z, column_symbol
    type(Operation_type) :: operation
    type(Atom_type) :: atom

    inquire(file=filename,exist=ex)
    if( .not. ex ) then
       write(*,'(a,a)') '# Error!: can not open file', trim(filename) 
       stop
    end if

    open(iunit, file=filename)

    do 
       read(iunit, '(a)', end=100) buf
       if( uncomment(buf) ) cycle
       read(buf,*) tag

       select case(tag)
       case("_cell_length_a")
          i = index(buf, "(")
          if( i == 0 ) then
             i = len(buf) + 1
          end if
          read(buf(1:i-1),*) tag, cif%cell%length_a

       case("_cell_length_b")
          i = index(buf, "(")
          if( i == 0 ) then
             i = len(buf) + 1
          end if
          read(buf(1:i-1),*) tag, cif%cell%length_b

       case("_cell_length_c")
          i = index(buf, "(")
          if( i == 0 ) then
             i = len(buf) + 1
          end if
          read(buf(1:i-1),*) tag, cif%cell%length_c

       case("_cell_angle_alpha")
          i = index(buf, "(")
          if( i == 0 ) then
             i = len(buf) + 1
          end if
          read(buf(1:i-1),*) tag, cif%cell%angle_alpha

       case("_cell_angle_beta")
          i = index(buf, "(")
          if( i == 0 ) then
             i = len(buf) + 1
          end if
          read(buf(1:i-1),*) tag, cif%cell%angle_beta

       case("_cell_angle_gamma")
          i = index(buf, "(")
          if( i == 0 ) then
             i = len(buf) + 1
          end if
          read(buf(1:i-1),*) tag, cif%cell%angle_gamma

       case("loop_")
          loop_atom = .false.
          loop_symm = .false.
          column = 0
          column_xyz   = -1
          column_x = -1
          column_y = -1
          column_z = -1
          column_symbol = -1

          do
             read(iunit, '(a)', end=100) buf
             if( uncomment(buf) ) cycle
             read(buf,*) tag
  
             if( tag=="" .or. tag(1:1) /= "_" ) exit

             column = column + 1
             select case(tag)
             case("_symmetry_equiv_pos_as_xyz")
                column_xyz = column
                loop_symm = .true.

             case("_atom_site_fract_x")
                column_x = column
                loop_atom = .true.

             case("_atom_site_fract_y")
                column_y = column
                loop_atom = .true.

             case("_atom_site_fract_z")
                column_z = column
                loop_atom = .true.

             case("_atom_site_type_symbol")
                column_symbol = column
                loop_atom = .true.
             end select
          end do
          column_max = column


          if( loop_symm ) then
             do
                if( buf=="" .or. buf == "loop_" .or. &
                     buf(1:1) == "_" .or. buf(1:1) == "#" ) exit

                call split(value,buf)

                operation = Operation_construct(value(column_xyz))

                call push_back_Operation( cif%voperation, operation )

                read(iunit, '(a)', end=100) buf
                if( uncomment(buf) ) exit
             end do
             cycle
          end if

          if( loop_atom ) then
             do
                if( buf=="" .or. buf == "loop_" .or. &
                     buf(1:1) == "_" .or. buf(1:1) == "#" ) exit

                call split(value,buf)

                read(value(column_x),*) atom%coords(1)
                read(value(column_y),*) atom%coords(2)
                read(value(column_z),*) atom%coords(3)
                read(value(column_symbol),*) atom%symbol
                call push_back_Atom( cif%vatom_cif, atom )

                read(iunit, '(a)', end=100) buf
                if( uncomment(buf) ) exit
             end do
             cycle
          end if
       end select
    end do
100 continue
    close(iunit)
  end subroutine CIF_readData


  !!----
  !! check cell length and angle.
  !! construct conventional unit vectors.
  !!----
  subroutine CIF_constructVector
    real(8) :: cosa
    real(8) :: cosb, sinb
    real(8) :: cosc, sinc
    real(8), parameter :: eps = 1.0d-5

    if( cif%cell%angle_alpha == 90.0d0 ) then
       cosa = 0.0d0
    else
       cosa = cos(M_PI/180.0*cif%cell%angle_alpha)
    end if
    if( cif%cell%angle_beta  == 90.0d0 ) then
       cosb = 0.0d0
       sinb = 1.0d0
    else
       cosb = cos(M_PI/180.0*cif%cell%angle_beta )
       sinb = sin(M_PI/180.0*cif%cell%angle_beta )
    end if
    if( cif%cell%angle_gamma == 90.0d0 ) then
       cosc = 0.0d0
       sinc = 1.0d0
    else
       cosc = cos(M_PI/180.0*cif%cell%angle_gamma)
       sinc = sin(M_PI/180.0*cif%cell%angle_gamma)
    end if

    cif%unit_conv(1,1) = cif%cell%length_a
    cif%unit_conv(2,1) = 0.0d0
    cif%unit_conv(3,1) = 0.0d0

    cif%unit_conv(1,2) = cif%cell%length_b*cosc
    cif%unit_conv(2,2) = cif%cell%length_b*sinc
    cif%unit_conv(3,2) = 0.0d0

    cif%unit_conv(1,3) = cif%cell%length_c*cosb
    cif%unit_conv(2,3) = cif%cell%length_c*(cosa-cosb*cosc)/sinc
    cif%unit_conv(3,3) = &
         sqrt( cif%cell%length_c**2 &
         - cif%unit_conv(1,3)**2 &
         - cif%unit_conv(2,3)**2 )

  end subroutine CIF_constructVector


  !!----
  !! calc position of atoms in conventional cell.
  !! calc position of atoms in primitive cell.
  !!----
  subroutine CIF_constructAtom
    integer :: ia, n
    type(Atom_type) :: atom

    !!  calc position of atoms in conventional cell
    do ia=1, size(cif%vatom_cif)
       do n=1, size(cif%voperation)

          atom%symbol = cif%vatom_cif(ia)%symbol
          atom%coords = Operation_operate( cif%voperation(n), cif%vatom_cif(ia)%coords )
          if( .not. Atom_find( cif%vatom_conv, atom ) ) then
             call push_back_Atom( cif%vatom_conv, atom )
          end if
       end do
    end do

  end subroutine CIF_constructAtom


  !!----
  !! calc symmetrize matrix for atom motion
  !!----
  subroutine CIF_constructSymmetrize
    integer :: ia, n, ib
    type(Atom_type) :: atom
    integer, allocatable :: multiplicity(:)

    allocate( cif%symmetrize(size(cif%vatom_conv),size(cif%vatom_conv)) )
    allocate( multiplicity(size(cif%vatom_conv)) )

    do ia=1, size(cif%vatom_conv)
       do ib=1, size(cif%vatom_conv)
          cif%symmetrize(ia,ib)%transM(:,:) = 0.0d0
          cif%symmetrize(ia,ib)%transV(:) = 0.0d0
       end do
    end do

    do ib=1, size(cif%vatom_conv)
       multiplicity(ib) = 0
    end do

    do ia=1, size(cif%vatom_conv)
       do n=1, size(cif%voperation)
          atom%symbol = cif%vatom_conv(ia)%symbol
          atom%coords = Operation_operate( cif%voperation(n), cif%vatom_conv(ia)%coords )
          do ib=1, size(cif%vatom_conv)
             if( Atom_match(cif%vatom_conv(ib),atom) ) then
                ! atom [a] is translated to atom [b] by operation [n]
                cif%symmetrize(ia,ib)%transM(:,:) &
                     = cif%symmetrize(ia,ib)%transM(:,:) &
                     + cif%voperation(n)%transM(:,:)
                multiplicity(ib) = multiplicity(ib) + 1
                exit
             end if
          end do
       end do
    end do

    do ib=1, size(cif%vatom_conv)
       if( multiplicity(ib) == 0 ) cycle

       do ia=1, size(cif%vatom_conv)
          cif%symmetrize(ia,ib)%transM(:,:) &
               = cif%symmetrize(ia,ib)%transM(:,:)/multiplicity(ib)
       end do
    end do

!!$    write(*,*) "DEBUG: restriction matrices"
!!$    do ia=1, size(cif%vatom_conv)
!!$       do ib=1, size(cif%vatom_conv)
!!$          if( any( cif%symmetrize(ia,ib)%transM(:,:) /= 0.0d0 ) ) then
!!$             write(*,"('a,b',2i3)") ia, ib
!!$             write(*,"(3f6.2)") cif%symmetrize(ia,ib)%transM(:,1)
!!$             write(*,"(3f6.2)") cif%symmetrize(ia,ib)%transM(:,2)
!!$             write(*,"(3f6.2)") cif%symmetrize(ia,ib)%transM(:,3)
!!$          end if
!!$       end do
!!$    end do

    deallocate( multiplicity )

  end subroutine CIF_constructSymmetrize


  !!----
  !! set CIF values to OptSW global variables.
  !!----
  subroutine CIF_setOptSW
    use paramlist
    integer :: ia
    logical :: symbol_found(103)

    QMD%uv(:,:) = cif%unit_conv(:,:)*(1.0d0/bohr)

    QMD%bv(:,3) = Position_outer_product( QMD%uv(:,1), QMD%uv(:,2) )
    QMD%bv(:,1) = Position_outer_product( QMD%uv(:,2), QMD%uv(:,3) )
    QMD%bv(:,2) = Position_outer_product( QMD%uv(:,3), QMD%uv(:,1) )
    QMD%omega = Position_inner_product( QMD%bv(:,3), QMD%uv(:,3) )
    QMD%omegai=1.d0/QMD%omega
    QMD%bv(:,:) = QMD%bv(:,:)*QMD%omegai

    QMD%natom = size(cif%vatom_conv)

    allocate( QMD%ra(3,QMD%natom), QMD%rr(3,QMD%natom) )
    allocate( QMD%rro(3,QMD%natom,2), QMD%iposfix(QMD%natom) )
    allocate( QMD%frc(3,QMD%natom), QMD%frco(3,QMD%natom,2), QMD%vrr(3,QMD%natom) )
    QMD%frc(:,:) = 0.d0
    QMD%vrr(:,:) = 0.d0

    allocate( QMD%zatm(QMD%natom), QMD%katm(QMD%natom) )
    allocate( QMD%mass(QMD%natom), QMD%mfac(QMD%natom) )

    QMD%mconv = (1.6605655d-27/9.109534d-31) ! Hartree atomic units

    symbol_found(:) = .false.
    do ia=1, size(cif%vatom_conv)
       QMD%zatm(ia) = getAtomicNumber(cif%vatom_conv(ia)%symbol)
       QMD%katm(ia) = 1  !!??
       QMD%iposfix(ia) = 1 !!??
       QMD%rr(:,ia) = cif%vatom_conv(ia)%coords(:)
       QMD%ra(:,ia) = matmul(QMD%uv(:,:),QMD%rr(:,ia))

       symbol_found( QMD%zatm(ia) ) = .true.
    end do

    QMD%nkatm = count( symbol_found(:) )
  end subroutine CIF_setOptSW

  logical function CIF_canSymmetrize()
    CIF_canSymmetrize = allocated( cif%symmetrize )
  end function CIF_canSymmetrize

  !!----
  !! symmetrize motion direction
  !!----
  subroutine CIF_symmetrizeDirection( drr )
    use paramlist

    real(8), intent(inout) :: drr(:,:) ! in relative coordinates
    real(8), allocatable :: drr_restrict(:,:)
    integer :: ia, ib

    if( size(cif%vatom_conv) /= QMD%natom ) then
       write(*,'(a,a)') '# Error!: array size mismatch in symmetrizeDirection'
       stop
    end if

    allocate( drr_restrict(3,size(cif%vatom_conv)) )

!!$    write(*,*) "DEBUG: symmetrize directions"
!!$    write(*,*) "  before symmetrize"
!!$    do ia=1, size(cif%vatom_conv)
!!$       write(*,"(i3,3e12.4)") ia, drr(:,ia)
!!$    end do

    ! calc restrict direction in internal coordinates
    do ib=1, size(cif%vatom_conv)
       drr_restrict(:,ib) = 0.0d0
       do ia=1, size(cif%vatom_conv)
          drr_restrict(:,ib) = drr_restrict(:,ib) &
               + matmul( cif%symmetrize(ia,ib)%transM(:,:), drr(:,ia) )
       end do
    end do

    ! overwrite direction by restrict direction
    do ia=1, size(cif%vatom_conv)
       drr(:,ia) = drr_restrict(:,ia)
    end do

!!$    write(*,*) "  after symmetrize"
!!$    do ia=1, size(cif%vatom_conv)
!!$       write(*,"(i3,3e12.4)") ia, drr(:,ia)
!!$    end do

    deallocate( drr_restrict )

  end subroutine CIF_symmetrizeDirection

  subroutine CIF_load( filename )
    character(len=*), intent(in) :: filename

    call CIF_readData(filename)
    call CIF_constructVector
    call CIF_constructAtom
    call CIF_constructSymmetrize
    call CIF_setOptSW
    
  end subroutine CIF_load

end module cif_module
