module EmDee

use, intrinsic :: iso_c_binding
use c_binding_extra

implicit none

integer, parameter, private :: ib = c_int
integer, parameter, private :: rb = c_double

integer(ib), parameter, private :: NONE  = 0, &
                                   LJ    = 1, &
                                   SF_LJ = 2

integer(ib), parameter, private :: ndiv = 2

integer(ib), parameter, private :: nbcells = 62

integer(ib), parameter, private :: nb(3,nbcells) = reshape( [ &
   1, 0, 0,    2, 0, 0,   -2, 1, 0,   -1, 1, 0,    0, 1, 0,    1, 1, 0,    2, 1, 0,   -2, 2, 0,  &
  -1, 2, 0,    0, 2, 0,    1, 2, 0,    2, 2, 0,   -2,-2, 1,   -1,-2, 1,    0,-2, 1,    1,-2, 1,  &
   2,-2, 1,   -2,-1, 1,   -1,-1, 1,    0,-1, 1,    1,-1, 1,    2,-1, 1,   -2, 0, 1,   -1, 0, 1,  &
   0, 0, 1,    1, 0, 1,    2, 0, 1,   -2, 1, 1,   -1, 1, 1,    0, 1, 1,    1, 1, 1,    2, 1, 1,  &
  -2, 2, 1,   -1, 2, 1,    0, 2, 1,    1, 2, 1,    2, 2, 1,   -2,-2, 2,   -1,-2, 2,    0,-2, 2,  &
   1,-2, 2,    2,-2, 2,   -2,-1, 2,   -1,-1, 2,    0,-1, 2,    1,-1, 2,    2,-1, 2,   -2, 0, 2,  &
  -1, 0, 2,    0, 0, 2,    1, 0, 2,    2, 0, 2,   -2, 1, 2,   -1, 1, 2,    0, 1, 2,    1, 1, 2,  &
   2, 1, 2,   -2, 2, 2,   -1, 2, 2,    0, 2, 2,    1, 2, 2,    2, 2, 2 ], [3,nbcells] )

integer(ib), private :: extra = 2000

type, bind(C) :: tCell
  integer(ib) :: neighbor(nbcells)
end type tCell

type, bind(C) :: tList
  integer(c_int) :: nitems
  integer(c_int) :: count
  type(c_ptr)    :: first
  type(c_ptr)    :: last
  type(c_ptr)    :: item
end type tList

type, bind(C) :: tPairType
  integer(ib) :: model
  real(rb)    :: p1
  real(rb)    :: p2
  real(rb)    :: p3
  real(rb)    :: p4
end type tPairType

type, bind(C) :: tEmDee

  integer(ib) :: builds    ! Number of neighbor-list builds

  type(tList) :: neighbor
  type(tList) :: excluded

  real(rb)    :: time

  integer(ib) :: natoms    ! Number of atoms
  integer(ib) :: nx3       ! Three times the number of atoms
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

  real(rb)    :: Energy
  real(rb)    :: Virial

end type tEmDee

private :: make_cells, distribute, find_pairs, compute

contains

!---------------------------------------------------------------------------------------------------

  subroutine reallocate_list( me, size )
    type(tList),    intent(inout) :: me
    integer(c_int), intent(in)    :: size

    integer(c_int) :: n
    integer(c_int), pointer :: old(:), new(:)

    call c_f_pointer( me%item, old, [me%nitems] )
    allocate( new(size) )
    n = min(me%nitems,size)
    new(1:n) = old(1:n)
    deallocate( old )
    me%item = c_loc(new)
    me%nitems = size

  end subroutine reallocate_list

!---------------------------------------------------------------------------------------------------

  subroutine make_cells( me, M )
    type(tEmDee), intent(inout) :: me
    integer(ib),  intent(in)    :: M

    integer(ib) :: MM, Mm1, k, ix, iy, iz
    type(tCell),  pointer :: cell(:)

    call c_f_pointer( me%cell, cell, [me%maxcells] )
    MM = M*M
    Mm1 = M - 1
    me%mcells = M
    me%ncells = M*MM
    if (me%ncells > me%maxcells) then
      deallocate( cell )
      allocate( cell(me%ncells) )
      me%maxcells = me%ncells
      me%cell = c_loc(cell)
    end if
    forall ( k = 1:nbcells, ix = 0:Mm1, iy = 0:Mm1, iz = 0:Mm1 )
      cell(1+ix+iy*M+iz*MM)%neighbor(k) = 1+pbc(ix+nb(1,k))+pbc(iy+nb(2,k))*M+pbc(iz+nb(3,k))*MM
    end forall
    nullify( cell )

    contains

      pure integer(ib) function pbc( x )
        integer(ib), intent(in) :: x
        if (x < 0) then
          pbc = x + M
        else if (x >= M) then
          pbc = x - M
        else
          pbc = x
        end if
      end function pbc

  end subroutine make_cells

!---------------------------------------------------------------------------------------------------

  subroutine distribute( M, N, R, atom, ncells, previous, natoms )
    integer(ib), intent(in)  :: M, N, ncells
    real(rb),    intent(in)  :: R(3,N)
    integer(ib), intent(out) :: atom(N), previous(ncells), natoms(ncells)
    integer(ib) :: i, j, icell, ntotal, MM
    integer(ib) :: next(N), head(ncells)

    head = 0
    natoms = 0
    MM = M*M
    do i = 1, N
      icell = 1 + int(M*R(1,i),ib) + M*int(M*R(2,i),ib) + MM*int(M*R(3,i),ib)
      next(i) = head(icell)
      head(icell) = i
      natoms(icell) = natoms(icell) + 1
    end do

    ntotal = 0
    do icell = 1, ncells
      previous(icell) = ntotal
      j = head(icell)
      do while (j /= 0)
        ntotal = ntotal + 1
        atom(ntotal) = j
        j = next(j)
      end do
    end do

  end subroutine distribute

!---------------------------------------------------------------------------------------------------

  subroutine find_pairs( me, L )
    type(tEmDee), intent(inout) :: me
    real(rb),     intent(in)    :: L

    integer(ib) :: i, j, k, m, nc, icell, jcell, maxpairs, npairs, nlocal, kprev, mprev, itype
    real(rb)    :: invL, invL2, xRc2
    real(rb)    :: Rs(3,me%natoms), Ri(3), Rij(3)
    integer(ib) :: natoms(me%ncells), previous(me%ncells), atom(me%natoms)

    type(tCell),     pointer :: cell(:)
    integer(ib),     pointer :: first(:), last(:), neighbor(:), type(:)
    integer(ib),     pointer :: xfirst(:), xlast(:), excluded(:)
    real(rb),        pointer :: R(:,:)
    type(tPairType), pointer :: pair(:,:)

    integer(ib), allocatable :: ilist(:)
    logical :: include

    call c_f_pointer( me%R, R, [3,me%natoms] )
    call c_f_pointer( me%cell, cell, [me%ncells] )
    call c_f_pointer( me%neighbor%first, first, [me%natoms] )
    call c_f_pointer( me%neighbor%last, last, [me%natoms] )
    call c_f_pointer( me%neighbor%item, neighbor, [me%neighbor%nitems] )
    call c_f_pointer( me%excluded%first, xfirst, [me%natoms] )
    call c_f_pointer( me%excluded%last, xlast, [me%natoms] )
    call c_f_pointer( me%excluded%item, excluded, [me%excluded%nitems] )
    call c_f_pointer( me%type, type, [me%natoms] )
    call c_f_pointer( me%pairType, pair, [me%ntypes,me%ntypes] )

    invL = 1.0_rb/L
    invL2 = invL*invL
    xRc2 = me%xRcSq*invL2

    ! Distribute atoms over cells and save scaled coordinates:
    Rs = R*invL
    Rs = Rs - floor(Rs)
    call distribute( me%mcells, me%natoms, Rs, atom, me%ncells, previous, natoms )

    ! Safely allocate local array:
    nc = maxval(natoms)
    maxpairs = (nc*((2*nbcells + 1)*nc - 1))/2

    ! Sweep all cells to search for neighbors:
    npairs = 0
    do icell = 1, me%ncells

      if (me%neighbor%nitems < npairs + maxpairs) then
        call reallocate_list( me%neighbor, npairs + maxpairs + extra )
        call c_f_pointer( me%neighbor%item, neighbor, [me%neighbor%nitems] )
      end if

      nlocal = natoms(icell)
      kprev = previous(icell)
      do k = 1, nlocal
        i = atom(kprev + k)
        first(i) = npairs + 1
        itype = type(i)
        Ri = Rs(:,i)
        ilist = excluded(xfirst(i):xlast(i))
        do m = k + 1, nlocal
          j = atom(kprev + m)
          call check_pair()
        end do
        do nc = 1, nbcells
          jcell = cell(icell)%neighbor(nc)
          mprev = previous(jcell)
          do m = 1, natoms(jcell)
            j = atom(mprev + m)
            call check_pair()
          end do
        end do
        last(i) = npairs
      end do

    end do
    me%neighbor%count = npairs
    nullify( R, cell, first, last, neighbor, type, pair )
    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine check_pair()
        if (pair(itype,type(j))%model /= NONE) then
          Rij = Ri - Rs(:,j)
          Rij = Rij - nint(Rij)
          if (sum(Rij*Rij) < xRc2) then
            if (i < j) then
              include = all(ilist /= j)
            else
              include = all(excluded(xfirst(j):xlast(j)) /= i)
            end if
            if (include) then
              npairs = npairs + 1
              neighbor(npairs) = j
            end if
          end if
        end if
      end subroutine check_pair
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine check_pair_alternative()
        if (pair(itype,type(j))%model /= NONE) then
          if (i < j) then
            include = all(ilist /= j)
          else
            include = all(excluded(xfirst(j):xlast(j)) /= i)
          end if
          if (include) then
            Rij = Ri - Rs(:,j)
            Rij = Rij - nint(Rij)
            if (sum(Rij*Rij) < xRc2) then
              npairs = npairs + 1
              neighbor(npairs) = j
            end if
          end if
        end if
      end subroutine check_pair_alternative
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine find_pairs

!---------------------------------------------------------------------------------------------------

  subroutine compute( me, L )
    type(tEmDee), intent(inout) :: me
    real(rb),     intent(in)    :: L

    integer  :: i, j, k, itype
    real(rb) :: invL, invL2, Rc2, Epot, Virial
    real(rb) :: r2, invR2, Rij(3), Ri(3), Fi(3), Fij(3), Eij, Wij
    real(rb) :: Rs(3,me%natoms), Fs(3,me%natoms)
    type(tPairType), pointer :: ij
    integer(ib),     pointer :: type(:), first(:), last(:), neighbor(:)
    real(rb),        pointer :: R(:,:), F(:,:)
    type(tPairType), pointer :: pairType(:,:)

    call c_f_pointer( me%type, type, [me%natoms] )
    call c_f_pointer( me%neighbor%first, first, [me%natoms] )
    call c_f_pointer( me%neighbor%last, last, [me%natoms] )
    call c_f_pointer( me%neighbor%item, neighbor, [me%neighbor%count] )
    call c_f_pointer( me%R, R, [3,me%natoms] )
    call c_f_pointer( me%F, F, [3,me%natoms] )
    call c_f_pointer( me%pairType, pairType, [me%ntypes,me%ntypes] )

    Epot = 0.0_8
    Virial = 0.0_8
    Fs = 0.0_8
    invL = 1.0_rb/L
    invL2 = invL*invL
    Rs = R*invL
    Rc2 = me%RcSq*invL2
    do i = 1, me%natoms
      itype = type(i)
      Ri = Rs(:,i)
      Fi = 0.0_rb
      do k = first(i), last(i)
        j = neighbor(k)
        Rij = Ri - Rs(:,j)
        Rij = Rij - nint(Rij)
        r2 = sum(Rij*Rij)
        if (r2 < Rc2) then
          invR2 = invL2/r2
          call compute_pair
          Epot = Epot + Eij
          Virial = Virial + Wij
          Fij = Wij*invR2*Rij
          Fi = Fi + Fij
          Fs(:,j) = Fs(:,j) - Fij
        end if
      end do
      Fs(:,i) = Fs(:,i) + Fi
    end do
    me%Energy = Epot
    me%Virial = Virial/3.0_rb
    F = L*Fs

    nullify( type, first, last, neighbor, R, F, pairType )
    contains
      include "pair_compute.f90"
  end subroutine compute

!---------------------------------------------------------------------------------------------------

  function match( a, b ) result( c )
    integer(ib), intent(in) :: a(:), b(:)
    logical                 :: c(size(a))
    ! Array "a" has unique entries in arbitrary order
    ! Array "b" has unique entries in increasing order
    ! Array "c" is such that c(i) is true iff a(i) is found in b
    integer :: n, i, value, mid, lb, ub
    integer :: p(size(b))
    n = size(b)
    p = b
    do i = 1, size(a)
      value = a(i)
      lb = 1
      ub = n
      mid = (lb + ub)/2
      do while( (p(mid) /= value) .and. (ub > lb) )
        if (value > p(mid)) then
          lb = mid + 1
        else
          ub = mid - 1
        end if
        mid = (lb + ub)/2
      end do
      c(i) = p(mid) == value
      if (c(i)) then
        n = n - 1
        p(mid:n) = p(mid+1:n+1)
      end if
    end do
  end function match

!---------------------------------------------------------------------------------------------------

  subroutine md_initialize( md, rc, skin, atoms, types, type_index, coords, forces ) bind(C)
    type(c_ptr), value :: md
    real(rb),    value :: rc, skin
    integer(ib), value :: atoms, types
    type(c_ptr), value :: type_index, coords, forces

    type(tEmDee),    pointer :: me
    integer(ib),     pointer :: type_ptr(:)
    type(tPairType), pointer :: pairType(:,:)

    call c_f_pointer( md, me )

    me%Rc = rc
    me%RcSq = rc*rc
    me%xRc = rc + skin
    me%xRcSq = me%xRc**2
    me%skinSq = skin*skin

    me%natoms = atoms
    me%nx3 = 3*atoms

    me%R = coords
    me%F = forces

    me%P  = malloc_real( me%nx3 )
    me%R0 = malloc_real( me%nx3, value = 0.0_rb )

    me%ntypes = types
    if (c_associated(type_index)) then
      call c_f_pointer( type_index, type_ptr, [atoms] )
      me%type = malloc_int( atoms, array = type_ptr )
      nullify( type_ptr )
    else
      me%type = malloc_int( atoms, value = 1 )
    end if

    allocate( pairType(types,types) )
    me%pairType = c_loc(pairType)
    pairType%model = NONE
    nullify(pairType)

    me%mcells = 0
    me%ncells = 0
    me%maxcells = 0
    me%builds = 0

    me%cell = malloc_int( 0 )

    me%neighbor%nitems = extra
    me%neighbor%count  = 0
    me%neighbor%item   = malloc_int( me%neighbor%nitems )
    me%neighbor%first  = malloc_int( atoms )
    me%neighbor%last   = malloc_int( atoms )

    me%excluded%nitems = extra
    me%excluded%count  = 0
    me%excluded%item   = malloc_int( me%excluded%nitems )
    me%excluded%first  = malloc_int( atoms, value = 1_ib )
    me%excluded%last   = malloc_int( atoms, value = 0_ib )

    me%time = 0.0_rb

    nullify( me )
  end subroutine md_initialize

!---------------------------------------------------------------------------------------------------

  subroutine md_set_pair( md, i, j, model ) bind(C)
    type(c_ptr),     value :: md
    integer(ib),     value :: i, j
    type(tPairType), value :: model
    type(tEmDee),    pointer :: me
    type(tPairType), pointer :: pairType(:,:)
    call c_f_pointer( md, me )
    call c_f_pointer( me%pairType, pairType, [me%ntypes,me%ntypes] )
    pairType(i,j) = model
    pairType(j,i) = model
    nullify( me, pairType )
  end subroutine md_set_pair

!---------------------------------------------------------------------------------------------------

  include "pair_setup.f90"

!---------------------------------------------------------------------------------------------------

  subroutine md_exclude_pair( md, i, j ) bind(C)
    type(c_ptr), value :: md
    integer(ib), value :: i, j

    integer :: n
    type(tEmDee), pointer :: me
    integer(ib),  pointer :: first(:), last(:), item(:)

    if (i /= j) then
      call c_f_pointer( md, me )
      call c_f_pointer( me%excluded%first, first, [me%natoms] )
      call c_f_pointer( me%excluded%last,  last,  [me%natoms] )
      n = me%excluded%count
      if (n == me%excluded%nitems) call reallocate_list( me%excluded, n + extra )
      call c_f_pointer( me%excluded%item, item, [me%excluded%nitems] )
      call add_item( min(i,j), max(i,j) )
!      call add_item( i, j )
!      call add_item( j, i )
      me%excluded%count = n
      nullify( me, first, last, item)
    end if

    contains

      subroutine add_item( i, j )
        integer, intent(in) :: i, j
        integer :: start, end
        start = first(i)
        end = last(i)
        if ((end < start).or.(j > item(end))) then
          item(end+2:n+1) = item(end+1:n)
          item(end+1) = j
        else
          do while (j > item(start))
            start = start + 1
          end do
          if (j == item(start)) return
          item(start+1:n+1) = item(start:n)
          item(start) = j
        end if
        first(i+1:) = first(i+1:) + 1
        last(i:) = last(i:) + 1
        n = n + 1
      end subroutine add_item

  end subroutine md_exclude_pair

!---------------------------------------------------------------------------------------------------

  subroutine md_compute_forces( md, L ) bind(C)
    type(c_ptr), value :: md
    real(rb),    value :: L

    integer  :: i, M
    real(rb) :: maximum, next, deltaSq
    real(rb) :: ti, tf
    type(tEmDee), pointer :: me
    real(rb),     pointer :: R(:,:), R0(:,:)

    call c_f_pointer( md, me )
    call c_f_pointer( me%R, R, [3,me%natoms] )
    call c_f_pointer( me%R0, R0, [3,me%natoms] )

    call cpu_time( ti )
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
      if (M /= me%mcells) call make_cells( me, M )
      call find_pairs( me, L )
      R0 = R
      me%builds = me%builds + 1
    end if

    call compute( me, L )
    
    call cpu_time( tf )
    me%time = me%time + (tf - ti)

    nullify( me, R, R0 )
  end subroutine md_compute_forces

!---------------------------------------------------------------------------------------------------

end module EmDee
