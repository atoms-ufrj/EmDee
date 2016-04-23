module mEmDee

use, intrinsic :: iso_c_binding

implicit none

integer, parameter, private :: ib = c_int
integer, parameter, private :: rb = c_double

integer, parameter, private :: NONE             = 0, &
                               LJ               = 1, &
                               SHIFTED_FORCE_LJ = 2

integer, parameter, private :: ndiv = 2

integer, parameter, private :: nbcells = 62

integer(ib), parameter, private :: nb(3,nbcells) = reshape( [ &
                               1, 0, 0,     2, 0, 0,    -2, 1, 0,    -1, 1, 0,     0, 1, 0,    &
                               1, 1, 0,     2, 1, 0,    -2, 2, 0,    -1, 2, 0,     0, 2, 0,    &
                               1, 2, 0,     2, 2, 0,    -2,-2, 1,    -1,-2, 1,     0,-2, 1,    &
                               1,-2, 1,     2,-2, 1,    -2,-1, 1,    -1,-1, 1,     0,-1, 1,    &
                               1,-1, 1,     2,-1, 1,    -2, 0, 1,    -1, 0, 1,     0, 0, 1,    &
                               1, 0, 1,     2, 0, 1,    -2, 1, 1,    -1, 1, 1,     0, 1, 1,    &
                               1, 1, 1,     2, 1, 1,    -2, 2, 1,    -1, 2, 1,     0, 2, 1,    &
                               1, 2, 1,     2, 2, 1,    -2,-2, 2,    -1,-2, 2,     0,-2, 2,    &
                               1,-2, 2,     2,-2, 2,    -2,-1, 2,    -1,-1, 2,     0,-1, 2,    &
                               1,-1, 2,     2,-1, 2,    -2, 0, 2,    -1, 0, 2,     0, 0, 2,    &
                               1, 0, 2,     2, 0, 2,    -2, 1, 2,    -1, 1, 2,     0, 1, 2,    &
                               1, 1, 2,     2, 1, 2,    -2, 2, 2,    -1, 2, 2,     0, 2, 2,    &
                               1, 2, 2,     2, 2, 2 ], [3,nbcells] )

type, bind(C) :: tCell
  integer(ib) :: neighbor(nbcells)
end type tCell

type, bind(C) :: tPairType
  integer(ib) :: model
  real(rb)    :: p1
  real(rb)    :: p2
  real(rb)    :: p3
  real(rb)    :: p4
end type tPairType

type, bind(C) :: tEmDee

  integer(ib) :: builds    ! Number of neighbor-list builds
  type(c_ptr) :: first     ! First neighbor of each atom
  type(c_ptr) :: last      ! Last neighbor of each atom 
  type(c_ptr) :: neighbor  ! List of neighbors

  real(rb)    :: neighbor_time
  real(rb)    :: pair_time

  integer(ib) :: natoms    ! Number of atoms
  integer(ib) :: nx3       ! Three times the number of atoms
  integer(ib) :: npairs    ! Number of neighbor pairs
  integer(ib) :: maxpairs  ! Maximum number of neighbor pairs
  integer(ib) :: mcells    ! Number of cells at each dimension
  integer(ib) :: ncells    ! Total number of cells
  integer(ib) :: maxcells  ! Maximum number of cells

  real(rb)    :: Rc        ! Cut-off distance
  real(rb)    :: RcSq      ! Cut-off distance squared
  real(rb)    :: xRc       ! Extended cutoff distance (including skin)
  real(rb)    :: xRcSq     ! Extended cutoff distance squared
  real(rb)    :: skinSq    ! Square of the neighbor list skin width

  type(c_ptr) :: cell

  type(c_ptr) :: type      ! Atom types
  type(c_ptr) :: R0        ! Atom positions at list building
  type(c_ptr) :: R         ! Pointer to dynamic atom positions
  type(c_ptr) :: P         ! Pointer to dynamic atom momenta
  type(c_ptr) :: F

  integer(ib) :: ntypes
  type(c_ptr) :: pairType
  type(c_ptr) :: invmass

  real(rb)    :: Energy
  real(rb)    :: Virial

end type tEmDee

interface

!  subroutine md_initialize( me, rc, skin, atoms, types, type_index, mass ) bind(C)
!    import :: c_ptr, ib, rb
!    type(c_ptr),    value :: me
!    real(rb), value :: rc, skin
!    integer(ib), value :: atoms, types
!    type(c_ptr),    value :: type_index, mass
!  end subroutine md_initialize

!  subroutine md_set_lj( me, i, j, sigma, epsilon ) bind(C)
!    import :: c_ptr, ib, rb
!    type(c_ptr),    value :: me
!    integer(ib), value :: i, j
!    real(rb), value :: sigma, epsilon
!  end subroutine md_set_lj

  subroutine md_set_shifted_force_lj( me, i, j, sigma, epsilon ) bind(C)
    import :: c_ptr, ib, rb
    type(c_ptr),    value :: me
    integer(ib), value :: i, j
    real(rb), value :: sigma, epsilon
  end subroutine md_set_shifted_force_lj

!  subroutine md_upload( me, coords, momenta ) bind(C)
!    import :: c_ptr
!    type(c_ptr), value :: me, coords, momenta
!  end subroutine md_upload

!  subroutine md_download( me, coords, momenta, forces ) bind(C)
!    import :: c_ptr
!    type(c_ptr), value :: me, coords, momenta, forces
!  end subroutine md_download

!  subroutine md_change_coordinates( me, a, b ) bind(C)
!    import :: c_ptr, rb
!    type(c_ptr),    value :: me
!    real(rb), value :: a, b
!  end subroutine md_change_coordinates

!  subroutine md_change_momenta( me, a, b ) bind(C)
!    import :: c_ptr, rb
!    type(c_ptr),    value :: me
!    real(rb), value :: a, b
!  end subroutine md_change_momenta

!  subroutine md_compute_forces( me, L ) bind(C)
!    import :: c_ptr, rb
!    type(c_ptr),    value :: me
!    real(rb), value :: L
!  end subroutine md_compute_forces

end interface

contains

!---------------------------------------------------------------------------------------------------

  subroutine make_cells( md, M ) bind(C)
    type(c_ptr), value :: md
    integer(ib), value :: M
    integer(ib) :: MM, Mm1, k, ix, iy, iz, jx, jy, jz, kx, ky, kz
    type(tEmDee), pointer :: me
    type(tCell),  pointer :: cell(:)
    call c_f_pointer( md, me )
    MM = M*M
    Mm1 = M - 1
    me%mcells = M
    me%ncells = M*MM
    call c_f_pointer( me%cell, cell, [me%maxcells] )
    if (me%ncells > me%maxcells) then
      deallocate( cell )
      allocate( cell(me%ncells) )
      me%maxcells = me%ncells
      me%cell = c_loc(cell)
    end if
    do k = 1, nbcells
      kx = nb(1,k)
      ky = nb(2,k)
      kz = nb(3,k)
      do iz = 0, Mm1
        jz = iz + kx
        call pbc( jz )
        do iy = 0, Mm1
          jy = iy + ky
          call pbc( jy )
          do ix = 0, Mm1
            jx = ix + kz
            call pbc( jx )
            cell(1 + ix + iy*M + iz*MM)%neighbor(k) = 1 + jx + jy*M + jz*MM
          end do
        end do
      end do
    end do
    nullify( cell )
    contains
      pure subroutine pbc( x )
        integer, intent(inout) :: x
        if (x < 0) then
          x = x + M
        else if (x >= M) then
          x = x - M
        end if
      end subroutine pbc
  end subroutine make_cells

!---------------------------------------------------------------------------------------------------

  subroutine find_pairs( md, L ) bind(C)
    type(c_ptr), value :: md
    real(rb),    value :: L

    integer(ib) :: M, MM, ntypes, i, j, k, icell, nmax, maxpairs, npairs, nlocal, ntotal
    real(rb)    :: invL, xRcSq, Ri(3), Rij(3)
    type(tEmDee),    pointer :: me
    type(tCell),     pointer :: cell(:)
    integer(ib),     pointer :: first(:), last(:), neighbor(:), type(:)
    real(rb),        pointer :: R(:,:)
    type(tPairType), pointer :: pairType(:)
    integer(ib), allocatable :: next(:), natoms(:), head(:), atom(:), iaux(:)
    real(rb),    allocatable :: SR(:,:)

    call c_f_pointer( md, me )
    call c_f_pointer( me%R, R, [3,me%natoms] )
    call c_f_pointer( me%cell, cell, [me%ncells] )
    call c_f_pointer( me%first, first, [me%natoms] )
    call c_f_pointer( me%last, last, [me%natoms] )
    call c_f_pointer( me%neighbor, neighbor, [me%maxpairs] )
    call c_f_pointer( me%type, type, [me%natoms] )
    call c_f_pointer( me%pairType, pairType, [me%ntypes*me%ntypes] )

    M = me%mcells
    MM = M*M
    ntypes = me%ntypes
    invL = 1.0_rb/L
    xRcSq = me%xRcSq*invL*invL

    ! Distribute atoms over cells and save scaled coordinates:
    allocate( SR(3,me%natoms), next(me%natoms), natoms(me%ncells), head(me%ncells) )
    SR = R*invL
    SR = SR - floor(SR)
    head = 0
    natoms = 0
    do i = 1, me%natoms
      icell = 1 + int(M*SR(1,i)) + M*int(M*SR(2,i)) + MM*int(M*SR(3,i))
      next(i) = head(icell)
      head(icell) = i
      natoms(icell) = natoms(icell) + 1
    end do

    ! Safely allocate local arrays:
    nmax = maxval(natoms)
    maxpairs = (nmax*((2*nbcells + 1)*nmax - 1))/2
    allocate( atom(nmax*(nbcells + 1)) )

    ! Sweep all cells to search for neighbors:
    npairs = 0
    do icell = 1, me%ncells
      nlocal = natoms(icell)
      if (nlocal /= 0) then

        if (me%maxpairs < npairs + maxpairs) then
          allocate( iaux(me%maxpairs) )
          iaux = neighbor
          deallocate( neighbor )
          me%maxpairs = npairs + maxpairs
          allocate( neighbor(me%maxpairs) )
          neighbor(1:size(iaux)) = iaux
          deallocate( iaux )
          me%neighbor = c_loc(neighbor)
        end if

        ! Build list of atoms in current cell and neighbor cells:
        ntotal = 0
        j = head(icell)
        do while (j /= 0)
          ntotal = ntotal + 1
          atom(ntotal) = j
          j = next(j)
        end do
        do k = 1, nbcells
          j = head(cell(icell)%neighbor(k))
          do while (j /= 0)
            ntotal = ntotal + 1
            atom(ntotal) = j
            j = next(j)
          end do
        end do

        ! Search for neighbors and add them to the list:
        do k = 1, nlocal
          i = atom(k)
          first(i) = npairs  ! CHANGE THIS LINE
          Ri = SR(:,i)
          do m = k+1, ntotal
            j = atom(m)
            if (pairType(ntypes*type(i) + type(j) + 1)%model /= NONE) then ! CHANGE THIS LINE
              Rij = Ri - SR(:,j)
              Rij = Rij - nint(Rij)
              if (sum(Rij*Rij) <= xRcSq) then
                npairs = npairs + 1
                neighbor(npairs) = j - 1   ! CHANGE THIS LINE
              end if
            end if
          end do
          last(i) = npairs - 1  ! CHANGE THIS LINE
        end do
      end if
    end do
    me%npairs = npairs

    nullify( me, R, cell, first, last, neighbor, type, pairType )

  end subroutine find_pairs

!---------------------------------------------------------------------------------------------------

  subroutine handle_neighbor_list( md, L ) bind(C)
    type(c_ptr), value :: md
    real(rb),    value :: L
    integer :: i, M
    real(rb) :: maximum, next, deltaSq
    real(rb),     pointer :: R(:,:), R0(:,:)
    type(tEmDee), pointer :: me

    call c_f_pointer( md, me )
    call c_f_pointer( me%R, R, [3,me%natoms] )
    call c_f_pointer( me%R0, R0, [3,me%natoms] )

    maximum = sum((R(:,1) - R0(:,1))**2)
    next = maximum
    do i = 2, me%natoms
      deltaSq = sum((R(:,i) - R0(:,i))**2)
      if (deltaSq > maximum) then
        next = maximum
        maximum = deltaSq
      end if
    end do

    if (maximum + 2.0_rb*sqrt(maximum*next) + next > me%skinSq) then
      M = floor(ndiv*L/me%xRc)
      if (M < 5) then
        write(0,'("ERROR: simulation box is too small.")')
        stop
      end if
      if (M /= me%mcells) call make_cells( md, M )
      call find_pairs( md, L )
      R0 = R
      me%builds = me%builds + 1
    end if

    nullify( me, R, R0 )

  end subroutine handle_neighbor_list

!---------------------------------------------------------------------------------------------------

  type(c_ptr) function malloc_int( n, value, array )
    integer(ib), intent(in)           :: n
    integer(ib), intent(in), optional :: value, array(:)
    integer(ib), pointer :: ptr(:)
    allocate( ptr(n) )
    malloc_int = c_loc(ptr)
    if (present(array)) then
      ptr = array
    else if (present(value)) then
      ptr = value
    end if
    nullify( ptr )
  end function malloc_int

!---------------------------------------------------------------------------------------------------

  type(c_ptr) function malloc_real( n, value, array )
    integer(ib), intent(in)           :: n
    real(rb),    intent(in), optional :: value, array(:)
    real(rb), pointer :: ptr(:)
    allocate( ptr(n) )
    malloc_real = c_loc(ptr)
    if (present(array)) then
      ptr = array
    else if (present(value)) then
      ptr = value
    end if
    nullify( ptr )
  end function malloc_real

!---------------------------------------------------------------------------------------------------

  subroutine copy_real( from, to, n )
    type(c_ptr), intent(in) :: from, to
    integer,     intent(in) :: n
    real(rb), pointer :: F(:), T(:)
    call c_f_pointer( from, F, [n] )
    call c_f_pointer( to, T, [n] )
    T = F
    nullify( T, F )
  end subroutine copy_real

!---------------------------------------------------------------------------------------------------

  subroutine md_initialize( md, rc, skin, atoms, types, type_index, mass ) bind(C)
    type(c_ptr), value :: md
    real(rb),    value :: rc, skin
    integer(ib), value :: atoms, types
    type(c_ptr), value :: type_index, mass

    type(tEmDee),    pointer :: me
    integer(ib),     pointer :: type_ptr(:)
    real(rb),        pointer :: mass_ptr(:)
    type(tPairType), pointer :: pair(:)

    call c_f_pointer( md, me )

    me%Rc = rc
    me%RcSq = rc*rc
    me%xRc = rc + skin
    me%xRcSq = (rc + skin)*(rc + skin)
    me%skinSq = skin*skin

    me%natoms = atoms
    me%nx3 = 3*atoms

    me%first = malloc_int( atoms )
    me%last = malloc_int( atoms )
    me%R = malloc_real( me%nx3 )
    me%P = malloc_real( me%nx3 )
    me%F = malloc_real( me%nx3 )
    me%R0 = malloc_real( me%nx3, 0.0_rb )

    me%ntypes = types
    call c_f_pointer( mass, mass_ptr, [types] )
    if (c_associated(type_index)) then
      call c_f_pointer( type_index, type_ptr, [atoms] )
      me%type = malloc_int( atoms, array = type_ptr - 1 )   ! CORRECT THIS LINE
      me%invmass = malloc_real( atoms, array = 1.0_rb/mass_ptr(type_ptr) )
      nullify( type_ptr )
    else
      me%type = malloc_int( atoms, value = 0 )   ! CORRECT THIS LINE
      me%invmass = malloc_real( atoms, value = 1.0_rb/mass_ptr(1) )
    end if
    nullify( mass_ptr )

    allocate( pair(types*types) )
    me%pairType = c_loc(pair)
    pair%model = NONE
    nullify(pair)

    me%mcells = 0
    me%ncells = 0
    me%maxcells = 0
    me%builds = 0

    me%cell = malloc_int( 0 )

    me%maxpairs = 1000000
    me%neighbor = malloc_int( me%maxpairs )

    me%neighbor_time = 0.0_rb
    me%pair_time = 0.0_rb

    nullify( me )
  end subroutine md_initialize

!---------------------------------------------------------------------------------------------------

  subroutine md_upload( md, coords, momenta ) bind(C)
    type(c_ptr), value :: md, coords, momenta
    type(tEmDee), pointer :: me
    call c_f_pointer( md, me )
    if (c_associated(coords)) call copy_real( coords, me%R, me%nx3 )
    if (c_associated(momenta)) call copy_real( momenta, me%P, me%nx3 )
  end subroutine md_upload

!---------------------------------------------------------------------------------------------------

  subroutine md_download( md, coords, momenta, forces ) bind(C)
    type(c_ptr), value :: md, coords, momenta, forces
    type(tEmDee), pointer :: me
    call c_f_pointer( md, me )
    if (c_associated(coords)) call copy_real( me%R, coords, me%nx3 )
    if (c_associated(momenta)) call copy_real( me%P, momenta, me%nx3 )
    if (c_associated(forces)) call copy_real( me%F, forces, me%nx3 )
  end subroutine md_download

!---------------------------------------------------------------------------------------------------

  subroutine md_set_lj( md, i, j, sigma, epsilon ) bind(C)
    type(c_ptr), value :: md
    integer(ib), value :: i, j
    real(rb),    value :: sigma, epsilon

    integer :: ij(2)
    type(tEmDee),    pointer :: me
    type(tPairType), pointer :: pairType(:)

    call c_f_pointer( md, me )
    call c_f_pointer( me%pairType, pairType, [me%ntypes**2] )
    ij = [i + me%ntypes*(j-1), j + me%ntypes*(i-1)]
    pairType(ij)%model = LJ
    pairType(ij)%p1 = sigma*sigma
    pairType(ij)%p2 = 4.0*epsilon
    nullify( me, pairType )

  end subroutine md_set_lj

!---------------------------------------------------------------------------------------------------

  subroutine md_change_coordinates( md, a, b ) bind(C)
    type(c_ptr), value :: md
    real(rb),    value :: a, b

    integer :: i
    type(tEmDee), pointer :: me
    real(rb),     pointer :: R(:,:), P(:,:), invmass(:)

    call c_f_pointer( md, me )
    call c_f_pointer( me%R, R, [3,me%natoms] )
    call c_f_pointer( me%P, P, [3,me%natoms] )
    call c_f_pointer( me%invmass, invmass, [me%natoms] )
    forall (i=1:me%natoms) R(:,i) = a*R(:,i) + b*invmass(i)*P(:,i)
    nullify( me, R, P, invmass )

  end subroutine md_change_coordinates

!---------------------------------------------------------------------------------------------------

  subroutine md_change_momenta( md, a, b ) bind(C)
    type(c_ptr), value :: md
    real(rb),    value :: a, b

    type(tEmDee), pointer :: me
    real(rb),     pointer :: P(:), F(:)

    call c_f_pointer( md, me )
    call c_f_pointer( me%P, P, [me%nx3] )
    call c_f_pointer( me%F, F, [me%nx3] )
    P = a*P + b*F
    nullify( me, P, F )

  end subroutine md_change_momenta

!---------------------------------------------------------------------------------------------------

  subroutine md_compute_forces( md, L ) bind(C)
    type(c_ptr), value :: md
    real(rb),    value :: L

    integer  :: i, j, k
    real(rb) :: InvL, InvL2, Epot, Virial, Rc2, r2, Rij(3), Ri(3), InvR2, InvR12, InvR6, Fi(3), Fij(3), Eij, Wij
    real(rb) :: ti, tm, tf
    integer(ib),  pointer :: type(:), first(:), last(:), neighbor(:)
    real(rb), pointer :: me_R(:,:), me_F(:,:)
    real(rb), allocatable :: R(:,:), F(:,:)
    type(tEmDee), pointer :: me

    call cpu_time( ti )
    call handle_neighbor_list( md, L )
    call cpu_time( tm )
    call c_f_pointer( md, me )
    call c_f_pointer( me%type, type, [me%natoms] )
    call c_f_pointer( me%first, first, [me%natoms] )
    call c_f_pointer( me%last, last, [me%natoms] )
    call c_f_pointer( me%neighbor, neighbor, [me%npairs] )

    allocate( F(3,me%natoms), R(3,me%natoms) )

    call c_f_pointer( me%R, me_R, [3,me%natoms] )

    Epot = 0.0_8
    Virial = 0.0_8
    F = 0.0_8
    InvL = 1.0_rb/L
    InvL2 = InvL*InvL
    R = me_R*invL
    Rc2 = me%RcSq*InvL2

    do i = 1, me%natoms
      Ri = R(:,i)
!      Fi = 0.0_rb
      do k = first(i)+1, last(i)+1
        j = neighbor(k)+1
        Rij = Ri - R(:,j)
        Rij = Rij - nint(Rij)
        r2 = sum(Rij*Rij)
        if (r2 < Rc2) then
          InvR2 = InvL2/r2
          InvR6 = InvR2*InvR2*InvR2
          InvR12 = InvR6*InvR6
          Eij = InvR12 - InvR6
          Wij = InvR12 + Eij

          Epot = Epot + Eij
          Virial = Virial + Wij
          Fij = L*Wij*InvR2*Rij

!          Fi = Fi + Fij
          F(:,i) = F(:,i) + Fij
          F(:,j) = F(:,j) - Fij
        end if
      end do
!      F(:,i) = F(:,i) + Fi
    end do
    me%Energy = 4.0_8*Epot
    me%Virial = 8.0_8*Virial
    F = 24.0_8*F

    call c_f_pointer( me%F, me_F, [3,me%natoms] )
    me_F = F

    call cpu_time( tf )

    me%neighbor_time = me%neighbor_time + (tm - ti)
    me%pair_time = me%pair_time + (tf - tm)

    nullify( me, type, me_R, me_F, first, last, neighbor )

  end subroutine md_compute_forces

!---------------------------------------------------------------------------------------------------

end module mEmDee
