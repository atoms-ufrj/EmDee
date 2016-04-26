! TO DO: Employ and debug subroutine distribute_atoms

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

integer(ib), private :: extraPairs = 2000

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

  real(rb)    :: time

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

contains

!---------------------------------------------------------------------------------------------------

  subroutine make_cells( me, M ) bind(C)
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

  subroutine find_pairs_and_compute( me, L )
    type(tEmDee), intent(inout) :: me
    real(rb),     intent(in)    :: L

    integer(ib) :: i, j, k, m, icell, jcell, nmax, maxpairs, npairs, nlocal, ntotal, itype, nprev
    real(rb)    :: invL, invL2, xRc2, Rc2
    real(rb)    :: Ri(3), Rij(3), Fi(3), Fij(3), Epot, Virial, r2, invR2, Eij, Wij
    integer(ib) :: natoms(me%ncells), previous(me%ncells), cell_atom(me%natoms)
    real(rb)    :: Rs(3,me%natoms), Fs(3,me%natoms)
    type(tPairType), pointer :: ij
    type(tCell),     pointer :: cell(:)
    integer(ib),     pointer :: first(:), last(:), neighbor(:), type(:)
    real(rb),        pointer :: R(:,:), F(:,:)
    type(tPairType), pointer :: pairType(:,:)
    integer(ib), allocatable :: atom(:)

    call c_f_pointer( me%R, R, [3,me%natoms] )
    call c_f_pointer( me%F, F, [3,me%natoms] )
    call c_f_pointer( me%cell, cell, [me%ncells] )
    call c_f_pointer( me%first, first, [me%natoms] )
    call c_f_pointer( me%last, last, [me%natoms] )
    call c_f_pointer( me%neighbor, neighbor, [me%maxpairs] )
    call c_f_pointer( me%type, type, [me%natoms] )
    call c_f_pointer( me%pairType, pairType, [me%ntypes,me%ntypes] )

    invL = 1.0_rb/L
    invL2 = invL*invL
    xRc2 = me%xRcSq*invL2
    Rc2 = me%RcSq*invL2

    ! Distribute atoms over cells and save scaled coordinates:
    Rs = R*invL
    call distribute( me%mcells, me%natoms, Rs - floor(Rs), cell_atom, me%ncells, previous, natoms )

    ! Safely allocate local array:
    nmax = maxval(natoms)
    maxpairs = (nmax*((2*nbcells + 1)*nmax - 1))/2
    allocate( atom(nmax*(nbcells + 1)) )

    ! Sweep all cells to search for neighbors:
    Fs = 0.0_rb
    Epot = 0.0_rb
    Virial = 0.0_rb
    npairs = 0
    do icell = 1, me%ncells
      nlocal = natoms(icell)
      if (nlocal /= 0) then

        if (me%maxpairs < npairs + maxpairs) then
          call realloc_int( me%neighbor, me%maxpairs, npairs + maxpairs + extraPairs )
          call c_f_pointer( me%neighbor, neighbor, [me%maxpairs] )
        end if

        ntotal = 0
        nprev = previous(icell)
        do j = 1, nlocal
          ntotal = ntotal + 1
          atom(ntotal) = cell_atom(nprev+j)
        end do
        do k = 1, nbcells
          jcell = cell(icell)%neighbor(k)
          nprev = previous(jcell)
          do j = 1, natoms(jcell)
            ntotal = ntotal + 1
            atom(ntotal) = cell_atom(nprev+j)
          end do
        end do

        ! Search for neighbors and add them to the list:
        do k = 1, nlocal
          i = atom(k)
          first(i) = npairs + 1
          itype = type(i)
          Ri = Rs(:,i)
          Fi = 0.0_rb
          do m = k+1, ntotal
            j = atom(m)
            ij => pairType(itype,type(j))
            if (ij%model /= NONE) then
              Rij = Ri - Rs(:,j)
              Rij = Rij - nint(Rij)
              r2 = sum(Rij*Rij)
              if (r2 < xRc2) then
                npairs = npairs + 1
                neighbor(npairs) = j
                call compute_pair
              end if
            end if
          end do
          Fs(:,i) = Fs(:,i) + Fi
          last(i) = npairs
        end do

      end if
    end do
    me%npairs = npairs
    me%Energy = Epot
    me%Virial = Virial/3.0_rb
    F = L*Fs
    nullify( R, F, cell, first, last, neighbor, type, pairType )
    contains
      include "pair_compute.f90"
  end subroutine find_pairs_and_compute

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
    call c_f_pointer( me%first, first, [me%natoms] )
    call c_f_pointer( me%last, last, [me%natoms] )
    call c_f_pointer( me%neighbor, neighbor, [me%npairs] )
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
        call compute_pair
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

  subroutine md_initialize( md, rc, skin, atoms, types, type_index, mass ) bind(C)
    type(c_ptr), value :: md
    real(rb),    value :: rc, skin
    integer(ib), value :: atoms, types
    type(c_ptr), value :: type_index, mass

    type(tEmDee),    pointer :: me
    integer(ib),     pointer :: type_ptr(:)
    real(rb),        pointer :: mass_ptr(:)
    type(tPairType), pointer :: pairType(:,:)

    call c_f_pointer( md, me )

    me%Rc = rc
    me%RcSq = rc*rc
    me%xRc = rc + skin
    me%xRcSq = me%xRc**2
    me%skinSq = skin*skin

    me%natoms = atoms
    me%nx3 = 3*atoms

    me%first = malloc_int( atoms )
    me%last  = malloc_int( atoms )
    me%R  = malloc_real( me%nx3 )
    me%P  = malloc_real( me%nx3 )
    me%F  = malloc_real( me%nx3 )
    me%R0 = malloc_real( me%nx3, value = 0.0_rb )

    me%ntypes = types
    call c_f_pointer( mass, mass_ptr, [types] )
    if (c_associated(type_index)) then
      call c_f_pointer( type_index, type_ptr, [atoms] )
      me%type = malloc_int( atoms, array = type_ptr )
      me%invmass = malloc_real( atoms, array = 1.0_rb/mass_ptr(type_ptr) )
      nullify( type_ptr )
    else
      me%type = malloc_int( atoms, value = 1 )
      me%invmass = malloc_real( atoms, value = 1.0_rb/mass_ptr(1) )
    end if
    nullify( mass_ptr )

    allocate( pairType(types,types) )
    me%pairType = c_loc(pairType)
    pairType%model = NONE
    nullify(pairType)

    me%mcells = 0
    me%ncells = 0
    me%maxcells = 0
    me%builds = 0

    me%cell = malloc_int( 0 )

    me%maxpairs = extraPairs
    me%neighbor = malloc_int( me%maxpairs )

    me%time = 0.0_rb

    nullify( me )
  end subroutine md_initialize

!---------------------------------------------------------------------------------------------------

  subroutine md_upload( md, coords, momenta ) bind(C)
    type(c_ptr), value :: md, coords, momenta
    type(tEmDee), pointer :: me
    call c_f_pointer( md, me )
    if (c_associated(coords))  call copy_real( coords, me%R, me%nx3 )
    if (c_associated(momenta)) call copy_real( momenta, me%P, me%nx3 )
  end subroutine md_upload

!---------------------------------------------------------------------------------------------------

  subroutine md_download( md, coords, momenta, forces ) bind(C)
    type(c_ptr), value :: md, coords, momenta, forces
    type(tEmDee), pointer :: me
    call c_f_pointer( md, me )
    if (c_associated(coords))  call copy_real( me%R, coords, me%nx3 )
    if (c_associated(momenta)) call copy_real( me%P, momenta, me%nx3 )
    if (c_associated(forces))  call copy_real( me%F, forces, me%nx3 )
  end subroutine md_download

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
      call find_pairs_and_compute( me, L )
      R0 = R
      me%builds = me%builds + 1
    else
      call compute( me, L )
    end if
    call cpu_time( tf )
    me%time = me%time + (tf - ti)

    nullify( me, R, R0 )
  end subroutine md_compute_forces

!---------------------------------------------------------------------------------------------------

end module EmDee