!   This file is part of EmDee.
!
!    EmDee is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    EmDee is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with EmDee. If not, see <http://www.gnu.org/licenses/>.
!
!    Author: Charlles R. A. Abreu (abreu@eq.ufrj.br)
!            Applied Thermodynamics and Molecular Simulation
!            Federal University of Rio de Janeiro, Brazil

module EmDee

use c_binding
use lists
use models
use bonded_structs

implicit none

real(rb), parameter, private :: zero  = 0.0_rb,                 &
                                one   = 1.0_rb,                 &
                                third = 0.33333333333333333_rb, &
                                pi    = 3.14159265358979324_rb

integer(ib), parameter, private :: extra = 2000
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

type, bind(C) :: tCell
  integer(ib) :: neighbor(nbcells)
end type tCell

type, bind(C) :: tEmDee

  integer(ib) :: builds      ! Number of neighbor-list builds

  type(tList) :: neighbor    ! List of neighbor atoms for pair interaction calculations
  type(tList) :: excluded    ! List of atom pairs excluded from the neighbor list

  real(rb)    :: time        ! Total time taken in force calculations

  integer(ib) :: natoms      ! Number of atoms in the system
  integer(ib) :: nx3         ! Three times the number of atoms
  integer(ib) :: mcells      ! Number of cells at each dimension
  integer(ib) :: ncells      ! Total number of cells
  integer(ib) :: maxcells    ! Maximum number of cells

  real(rb)    :: Rc          ! Cut-off distance
  real(rb)    :: RcSq        ! Cut-off distance squared
  real(rb)    :: xRc         ! Extended cutoff distance (including skin)
  real(rb)    :: xRcSq       ! Extended cutoff distance squared
  real(rb)    :: skinSq      ! Square of the neighbor list skin width

  type(c_ptr) :: cell        ! Array containing all neighbor cells of each cell

  type(c_ptr) :: type        ! The type of each atom
  type(c_ptr) :: R0          ! The position of each atom at the latest neighbor list building
  type(c_ptr) :: R           ! Pointer to the coordinate of each atom
  type(c_ptr) :: F           ! Pointer to the resultant force over each atom
  type(c_ptr) :: charge      ! Pointer to the electric charge of each atom

  integer(ib) :: ntypes      ! Number of atom types
  type(c_ptr) :: pairType    ! Model and parameters of each type of atom pair
  real(rb)    :: pairEnergy  ! Potential energy due to pair interactions
  real(rb)    :: pairVirial  ! Internal virial due to pair interactions

  integer(ib) :: nbonds      ! Number of chemical bonds in the system
  integer(ib) :: maxbonds    ! Maximum number of stored chemical bonds
  type(c_ptr) :: bond        ! Pointer to the list of chemical bonds
  real(rb)    :: bondEnergy  ! Potential energy due to bond stretching
  real(rb)    :: bondVirial  ! Internal virial due to bond stretching

  integer(ib) :: nangles     ! Number of bond angles in the system
  integer(ib) :: maxangles   ! Maximum number of stored bond angles
  type(c_ptr) :: angle       ! Pointer to the list of bond angles
  real(rb)    :: angleEnergy ! Potential energy due to angle bending
  real(rb)    :: angleVirial ! Internal virial due to angle bending

  real(rb)    :: Energy      ! Total potential energy of the system
  real(rb)    :: Virial      ! Total internal virial of the system

end type tEmDee

private :: reallocate_list, add_bonded_struc, make_cells, distribute_atoms, &
           compute_bonds, compute_angles, find_pairs_and_compute, compute_pairs
           

contains

!===================================================================================================
!                                L I B R A R Y   P R O C E D U R E S
!===================================================================================================

  subroutine md_initialize( md, rc, skin, atoms, types, indices, charges, coords, forces ) bind(C)
    type(c_ptr), value :: md
    real(rb),    value :: rc, skin
    integer(ib), value :: atoms, types
    type(c_ptr), value :: indices, charges, coords, forces

    type(tEmDee),    pointer :: me
    integer(ib),     pointer :: type_ptr(:)
    real(rb),        pointer :: charge_ptr(:)
    type(tModelPtr), pointer :: pairType(:,:)

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
    me%ntypes = types
    me%mcells = 0
    me%ncells = 0
    me%maxcells = 0
    me%nbonds = 0
    me%maxbonds = 0
    me%bondEnergy = zero
    me%bondVirial = zero
    me%nangles = 0
    me%maxangles = 0
    me%angleEnergy = zero
    me%angleVirial = zero
    me%builds = 0
    me%time = zero
    me%R0 = malloc_real( me%nx3, value = zero )
    me%cell = malloc_int( 0 )
    me%bond = malloc_int( 0 )
    me%angle = malloc_int( 0 )

    call allocate_list( me%neighbor, extra, atoms )
    call allocate_list( me%excluded, extra, atoms )

    if (c_associated(indices)) then
      call c_f_pointer( indices, type_ptr, [atoms] )
      me%type = malloc_int( atoms, array = type_ptr )
    else
      me%type = malloc_int( atoms, value = 1 )
    end if

    if (c_associated(charges)) then
      call c_f_pointer( charges, charge_ptr, [atoms] )
      me%charge = malloc_real( atoms, array = charge_ptr )
    else
      me%charge = malloc_real( atoms, value = zero )
    end if

    allocate( pairType(types,types) )
    me%pairType = c_loc(pairType(1,1))

    nullify( me, type_ptr, charge_ptr, pairType )
  end subroutine md_initialize

!---------------------------------------------------------------------------------------------------

  subroutine md_set_pair( md, itype, jtype, model ) bind(C)
    type(c_ptr), value :: md
    integer(ib), value :: itype, jtype
    type(c_ptr), value :: model

    type(tEmDee),    pointer :: me
    type(tModelPtr), pointer :: pairType(:,:)

    call c_f_pointer( md, me )
    call c_f_pointer( me%pairType, pairType, [me%ntypes,me%ntypes] )
    call c_f_pointer( model, pairType(itype,jtype)%model )
    call c_f_pointer( model, pairType(jtype,itype)%model )

    nullify( me, pairType )
  end subroutine md_set_pair

!---------------------------------------------------------------------------------------------------

  subroutine md_add_bond( md, i, j, model ) bind(C)
    type(c_ptr), value :: md
    integer(ib), value :: i, j
    type(c_ptr), value :: model

    type(tEmDee), pointer :: me

    call c_f_pointer( md, me )
    call add_bonded_struc( me%bond, me%nbonds, me%maxbonds, i, j, 0, 0, model )
    call md_exclude_pair( md, i, j )

    nullify( me )
  end subroutine md_add_bond

!---------------------------------------------------------------------------------------------------

  subroutine md_add_angle( md, i, j, k, model ) bind(C)
    type(c_ptr), value :: md
    integer(ib), value :: i, j, k
    type(c_ptr), value :: model

    type(tEmDee), pointer :: me

    call c_f_pointer( md, me )
    call add_bonded_struc( me%angle, me%nangles, me%maxangles, i, j, k, 0, model )
    call md_exclude_pair( md, i, j )
    call md_exclude_pair( md, j, k )

    nullify( me )
  end subroutine md_add_angle

!---------------------------------------------------------------------------------------------------

  subroutine md_exclude_pair( md, i, j ) bind(C)
    type(c_ptr), value :: md
    integer(ib), value :: i, j

    integer(ib) :: n
    type(tEmDee), pointer :: me
    integer(ib),  pointer :: first(:), last(:), item(:)

    call c_f_pointer( md, me )
    if ((i > 0).and.(i <= me%natoms).and.(j > 0).and.(j <= me%natoms).and.(i /= j)) then
      call c_f_pointer( me%excluded%first, first, [me%natoms] )
      call c_f_pointer( me%excluded%last,  last,  [me%natoms] )
      n = me%excluded%count
      if (n == me%excluded%nitems) call reallocate_list( me%excluded, n + extra )
      call c_f_pointer( me%excluded%item, item, [me%excluded%nitems] )
      call add_item( min(i,j), max(i,j) )
      me%excluded%count = n
      nullify( first, last, item)
    end if

    nullify( me )
    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine md_exclude_pair

!---------------------------------------------------------------------------------------------------

  subroutine md_compute_forces( md, L ) bind(C)
    type(c_ptr), value :: md
    real(rb),    value :: L

    integer(ib) :: M
    real(rb)    :: time
    type(tEmDee), pointer :: me

    call c_f_pointer( md, me )

    call cpu_time( time )
    me%time = me%time - time

    if (maximum_approach_sq( me ) > me%skinSq) then
      M = floor(ndiv*L/me%xRc)
      if (M < 5) then
        write(0,'("ERROR: simulation box exceedingly small.")')
        stop
      end if
      if (M /= me%mcells) call make_cells( me, M )
      call find_pairs_and_compute( me, L )
      call copy_real( me%R, me%R0, me%nx3 )
      me%builds = me%builds + 1
    else
      call compute_pairs( me, L )
    end if

    if (me%nbonds  /= 0) call compute_bonds( me, L )
    if (me%nangles /= 0) call compute_angles( me, L )

    me%Energy = me%pairEnergy + me%bondEnergy + me%angleEnergy
    me%Virial = me%pairVirial + me%bondVirial + me%angleVirial

    call cpu_time( time )
    me%time = me%time + time

    nullify( me )
  end subroutine md_compute_forces

!===================================================================================================
!                              A U X I L I A R Y   P R O C E D U R E S
!===================================================================================================

  real(rb) function maximum_approach_sq( me )
    type(tEmDee), intent(in) :: me

    integer(ib) :: i
    real(rb)    :: maximum, next, deltaSq
    real(rb), pointer :: R(:,:), R0(:,:)

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
    maximum_approach_sq = maximum + 2*sqrt(maximum*next) + next

    nullify( R, R0 )
  end function maximum_approach_sq

!---------------------------------------------------------------------------------------------------

  subroutine make_cells( me, M )
    type(tEmDee), intent(inout) :: me
    integer(ib),  intent(in)    :: M

    integer(ib) :: MM, Mm1, k, ix, iy, iz
    type(tCell), pointer :: cell(:)

    call c_f_pointer( me%cell, cell, [me%maxcells] )
    MM = M*M
    Mm1 = M - 1
    me%mcells = M
    me%ncells = M*MM
    if (me%ncells > me%maxcells) then
      deallocate( cell )
      allocate( cell(me%ncells) )
      me%maxcells = me%ncells
      me%cell = c_loc(cell(1))
    end if
    forall ( k = 1:nbcells, ix = 0:Mm1, iy = 0:Mm1, iz = 0:Mm1 )
      cell(1+ix+iy*M+iz*MM)%neighbor(k) = 1+pbc(ix+nb(1,k))+pbc(iy+nb(2,k))*M+pbc(iz+nb(3,k))*MM
    end forall
    nullify( cell )

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine make_cells

!---------------------------------------------------------------------------------------------------

  subroutine distribute_atoms( M, N, R, atom, ncells, previous, natoms )
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

  end subroutine distribute_atoms

!---------------------------------------------------------------------------------------------------

  subroutine compute_bonds( me, L )
    type(tEmDee), intent(inout) :: me
    real(rb),     intent(in)    :: L

    integer(ib) :: i, j, m
    real(rb)    :: invL, d, E, mdEdr, Energy, Virial
    real(rb)    :: Rij(3), Fij(3)

    type(tModel),  pointer :: model
    real(rb),      pointer :: R(:,:), F(:,:)
    type(tStruct), pointer :: bond(:)

    call c_f_pointer( me%bond, bond, [me%nbonds] )
    call c_f_pointer( me%R, R, [3,me%natoms] )
    call c_f_pointer( me%F, F, [3,me%natoms] )

    invL = one/L
    Energy = zero
    Virial = zero
    do m = 1, me%nbonds
      i = bond(m)%i
      j = bond(m)%j
      Rij = R(:,i) - R(:,j)
      Rij = Rij - L*nint(invL*Rij)
      d = sqrt(sum(Rij*Rij))
      model => bond(m)%model
      call compute_bond
      Energy = Energy + E
      Virial = Virial + mdEdr*d
      Fij = mdEdr*Rij/d
      F(:,i) = F(:,i) + Fij
      F(:,j) = F(:,j) - Fij
    end do
    me%bondEnergy = Energy
    me%bondVirial = third*Virial

    nullify( bond, R, F )
    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      include "compute_bond.f90"
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine compute_bonds

!---------------------------------------------------------------------------------------------------

  subroutine compute_angles( me, L )
    type(tEmDee), intent(inout) :: me
    real(rb),     intent(in)    :: L

    integer(ib) :: i, j, k, m
    real(rb)    :: invL, Energy, Virial, RijSq, RkjSq, dot, cross, theta, invDot, E, dEdtheta, Fa
    real(rb)    :: Rj(3), Rij(3), Rkj(3), Fi(3), Fk(3)

    type(tModel),  pointer :: model
    real(rb),      pointer :: R(:,:), F(:,:)
    type(tStruct), pointer :: angle(:)

    call c_f_pointer( me%angle, angle, [me%nangles] )
    call c_f_pointer( me%R, R, [3,me%natoms] )
    call c_f_pointer( me%F, F, [3,me%natoms] )

    invL = one/L
    Energy = zero
    Virial = zero
    do m = 1, me%nangles
      i = angle(m)%i
      j = angle(m)%j
      k = angle(m)%k
      Rj = R(:,j)
      Rij = R(:,i) - Rj
      Rkj = R(:,k) - Rj
      Rij = Rij - L*nint(invL*Rij)
      Rkj = Rkj - L*nint(invL*Rkj)
      RijSq = sum(Rij*Rij)
      RkjSq = sum(Rkj*Rkj)
      dot = sum(Rij*Rkj)
      cross = sqrt(RijSq*RkjSq - dot*dot)
      theta = pi + atan2(cross,dot)
      model => angle(m)%model
      call compute_angle()
      Fa = dEdtheta*dot/cross
      Fi = Fa*(Rkj*invDot - Rij/RijSq)
      Fk = Fa*(Rij*invDot - Rkj/RkjSq)
      F(:,i) = F(:,i) + Fi
      F(:,k) = F(:,k) + Fk
      F(:,j) = F(:,j) - (Fi + Fk)
      Energy = Energy + E
      Virial = Virial + sum(Fi*Rij + Fk*Rkj)
    end do
    me%angleEnergy = Energy
    me%angleVirial = third*Virial
    nullify( R, F, angle )
    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      include "compute_angle.f90"
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine compute_angles

!---------------------------------------------------------------------------------------------------

  subroutine find_pairs_and_compute( me, L )
    type(tEmDee), intent(inout) :: me
    real(rb),     intent(in)    :: L

    integer(ib) :: i, j, k, m, n, icell, jcell, maxpairs, npairs, nlocal, kprev, mprev, itype
    real(rb)    :: invL, invL2, xRc2, Rc2, r2, invR2, Epot, Virial, Eij, Wij, icharge
    logical     :: include
    integer(ib) :: natoms(me%ncells), previous(me%ncells), atom(me%natoms)
    real(rb)    :: Rs(3,me%natoms), Ri(3), Rij(3), Fi(3), Fij(3)

    integer(ib), allocatable :: ilist(:)

    type(tModel),    pointer :: model
    type(tCell),     pointer :: cell(:)
    integer(ib),     pointer :: first(:), last(:), neighbor(:), type(:)
    integer(ib),     pointer :: xfirst(:), xlast(:), excluded(:)
    real(rb),        pointer :: R(:,:), F(:,:), charge(:)
    type(tModelPtr), pointer :: pairType(:,:)

    call c_f_pointer( me%cell, cell, [me%ncells] )
    call c_f_pointer( me%neighbor%first, first, [me%natoms] )
    call c_f_pointer( me%neighbor%last, last, [me%natoms] )
    call c_f_pointer( me%neighbor%item, neighbor, [me%neighbor%nitems] )
    call c_f_pointer( me%excluded%first, xfirst, [me%natoms] )
    call c_f_pointer( me%excluded%last, xlast, [me%natoms] )
    call c_f_pointer( me%excluded%item, excluded, [me%excluded%nitems] )
    call c_f_pointer( me%type, type, [me%natoms] )
    call c_f_pointer( me%R, R, [3,me%natoms] )
    call c_f_pointer( me%F, F, [3,me%natoms] )
    call c_f_pointer( me%charge, charge, [me%natoms] )
    call c_f_pointer( me%pairType, pairType, [me%ntypes,me%ntypes] )

    invL = one/L
    invL2 = invL*invL
    xRc2 = me%xRcSq*invL2
    Rc2 = me%RcSq*invL2
    Rs = R*invL
    Rs = Rs - floor(Rs)

    call distribute_atoms( me%mcells, me%natoms, Rs, atom, me%ncells, previous, natoms )

    n = maxval(natoms)
    maxpairs = (n*((2*nbcells + 1)*n - 1))/2

    Epot = zero
    Virial = zero
    F = zero
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
        icharge = charge(i)
        Ri = Rs(:,i)
        Fi = zero
        ilist = excluded(xfirst(i):xlast(i))
        do m = k + 1, nlocal
          j = atom(kprev + m)
          call check_pair()
        end do
        do n = 1, nbcells
          jcell = cell(icell)%neighbor(n)
          mprev = previous(jcell)
          do m = 1, natoms(jcell)
            j = atom(mprev + m)
            call check_pair()
          end do
        end do
        F(:,i) = F(:,i) + Fi
        last(i) = npairs
      end do

    end do
    me%neighbor%count = npairs
    me%pairEnergy = Epot
    me%pairVirial = third*Virial
    F = L*F
    nullify( cell, first, last, neighbor, type, R, F, charge, pairType )
    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine check_pair()
        model => pairType(itype,type(j))%model
        if (associated(model)) then
          Rij = Ri - Rs(:,j)
          Rij = Rij - nint(Rij)
          r2 = sum(Rij*Rij)
          if (r2 < xRc2) then
            if (i < j) then
              include = all(ilist /= j)
            else
              include = all(excluded(xfirst(j):xlast(j)) /= i)
            end if
            if (include) then
              npairs = npairs + 1
              neighbor(npairs) = j
              if (r2 < Rc2) then
                invR2 = invL2/r2
                call compute_pair()
                Epot = Epot + Eij
                Virial = Virial + Wij
                Fij = Wij*invR2*Rij
                Fi = Fi + Fij
                F(:,j) = F(:,j) - Fij
              end if
            end if
          end if
        end if
      end subroutine check_pair
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      include "compute_pair.f90"
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine find_pairs_and_compute

!---------------------------------------------------------------------------------------------------

  subroutine compute_pairs( me, L )
    type(tEmDee), intent(inout) :: me
    real(rb),     intent(in)    :: L

    integer  :: i, j, k, itype
    real(rb) :: invL, invL2, Rc2, r2, invR2, Epot, Virial, Eij, Wij, icharge
    real(rb) :: Rs(3,me%natoms), Rij(3), Ri(3), Fi(3), Fij(3)

    type(tModel),    pointer :: model
    integer(ib),     pointer :: type(:), first(:), last(:), neighbor(:)
    real(rb),        pointer :: R(:,:), F(:,:), charge(:)
    type(tModelPtr), pointer :: pairType(:,:)

    call c_f_pointer( me%type, type, [me%natoms] )
    call c_f_pointer( me%R, R, [3,me%natoms] )
    call c_f_pointer( me%F, F, [3,me%natoms] )
    call c_f_pointer( me%charge, charge, [me%natoms] )
    call c_f_pointer( me%neighbor%first, first, [me%natoms] )
    call c_f_pointer( me%neighbor%last, last, [me%natoms] )
    call c_f_pointer( me%neighbor%item, neighbor, [me%neighbor%count] )
    call c_f_pointer( me%pairType, pairType, [me%ntypes,me%ntypes] )

    invL = one/L
    invL2 = invL*invL
    Rc2 = me%RcSq*invL2
    Rs = R*invL
    Epot = zero
    Virial = zero
    F = zero
    do i = 1, me%natoms
      itype = type(i)
      Ri = Rs(:,i)
      Fi = zero
      icharge = charge(i)
      do k = first(i), last(i)
        j = neighbor(k)
        Rij = Ri - Rs(:,j)
        Rij = Rij - nint(Rij)
        r2 = sum(Rij*Rij)
        if (r2 < Rc2) then
          invR2 = invL2/r2
          model => pairType(itype,type(j))%model
          call compute_pair()
          Epot = Epot + Eij
          Virial = Virial + Wij
          Fij = Wij*invR2*Rij
          Fi = Fi + Fij
          F(:,j) = F(:,j) - Fij
        end if
      end do
      F(:,i) = F(:,i) + Fi
    end do
    me%pairEnergy = Epot
    me%pairVirial = third*Virial
    F = L*F

    nullify( type, R, F, charge, first, last, neighbor, pairType )
    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      include "compute_pair.f90"
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine compute_pairs

!---------------------------------------------------------------------------------------------------

end module EmDee
