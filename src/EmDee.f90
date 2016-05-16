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
                                pi    = 3.14159265358979324_rb, &
                                piBy2 = 0.5_rb*pi

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

  type(tStructData) :: bond
  type(tStructData) :: angle
  type(tStructData) :: dihedral

  real(rb)    :: Energy      ! Total potential energy of the system
  real(rb)    :: Virial      ! Total internal virial of the system

end type tEmDee

private :: reallocate_list, make_cells, distribute_atoms, find_pairs_and_compute, &
           compute_pairs, compute_bonds, compute_angles, compute_dihedrals

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
    me%bond%number = 0
    me%bond%max = 0
    me%bond%energy = zero
    me%bond%virial = zero
    me%angle%number = 0
    me%angle%max = 0
    me%angle%energy = zero
    me%angle%virial = zero
    me%builds = 0
    me%time = zero
    me%R0 = malloc_real( me%nx3, value = zero )
    me%cell = malloc_int( 0 )
    me%bond%list = malloc_int( 0 )
    me%angle%list = malloc_int( 0 )

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
    call add_bonded_struc( me%bond, i, j, 0, 0, model )
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
    call add_bonded_struc( me%angle, i, j, k, 0, model )
    call md_exclude_pair( md, i, j )
    call md_exclude_pair( md, i, k )
    call md_exclude_pair( md, j, k )

    nullify( me )
  end subroutine md_add_angle

!---------------------------------------------------------------------------------------------------

  subroutine md_add_dihedral( md, i, j, k, l, model ) bind(C)
    type(c_ptr), value :: md
    integer(ib), value :: i, j, k, l
    type(c_ptr), value :: model

    type(tEmDee), pointer :: me

    call c_f_pointer( md, me )
    call add_bonded_struc( me%dihedral, i, j, k, l, model )
    call md_exclude_pair( md, i, j )
    call md_exclude_pair( md, i, k )
    call md_exclude_pair( md, i, l )
    call md_exclude_pair( md, j, k )
    call md_exclude_pair( md, j, l )
    call md_exclude_pair( md, k, l )

    nullify( me )
  end subroutine md_add_dihedral

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

    if (me%bond%number /= 0) call compute_bonds( me, L )
    if (me%angle%number /= 0) call compute_angles( me, L )
    if (me%dihedral%number /= 0) call compute_dihedrals( me, L )

    me%Energy = me%pairEnergy + me%bond%energy + me%angle%energy
    me%Virial = me%pairVirial + me%bond%virial + me%angle%virial

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

    call c_f_pointer( me%bond%list, bond, [me%bond%number] )
    call c_f_pointer( me%R, R, [3,me%natoms] )
    call c_f_pointer( me%F, F, [3,me%natoms] )

    invL = one/L
    Energy = zero
    Virial = zero
    do m = 1, me%bond%number
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
    me%bond%energy = Energy
    me%bond%virial = third*Virial

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
    real(rb)    :: invL, Energy, Virial, aa, bb, ab, axb, theta, Ea, Fa
    real(rb)    :: Rj(3), Fi(3), Fk(3), a(3), b(3)

    type(tModel),  pointer :: model
    real(rb),      pointer :: R(:,:), F(:,:)
    type(tStruct), pointer :: angle(:)

    call c_f_pointer( me%angle%list, angle, [me%angle%number] )
    call c_f_pointer( me%R, R, [3,me%natoms] )
    call c_f_pointer( me%F, F, [3,me%natoms] )

    invL = one/L
    Energy = zero
    Virial = zero
    do m = 1, me%angle%number
      i = angle(m)%i
      j = angle(m)%j
      k = angle(m)%k
      model => angle(m)%model
      Rj = R(:,j)
      a = R(:,i) - Rj
      b = R(:,k) - Rj
      a = a - L*nint(invL*a)
      b = b - L*nint(invL*b)
      aa = sum(a*a)
      bb = sum(b*b)
      ab = sum(a*b)
      axb = sqrt(aa*bb - ab*ab)
      theta = atan2(axb,ab)
      call compute_angle()
      Fa = Fa/axb
      Fi = Fa*(b - (ab/aa)*a)
      Fk = Fa*(a - (ab/bb)*b)
      F(:,i) = F(:,i) + Fi
      F(:,k) = F(:,k) + Fk
      F(:,j) = F(:,j) - (Fi + Fk)
      Energy = Energy + Ea
      Virial = Virial + sum(Fi*a + Fk*b)
    end do
    me%angle%energy = Energy
    me%angle%virial = third*Virial
    nullify( R, F, angle )
    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      include "compute_angle.f90"
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine compute_angles

!---------------------------------------------------------------------------------------------------

  subroutine compute_dihedrals( me, L )
    type(tEmDee), intent(inout) :: me
    real(rb),     intent(in)    :: L

    integer(ib) :: m
    real(rb)    :: invL, Energy, Virial, Ed, Fd, r2, invR2, Eij, Wij, icharge, jcharge
    real(rb)    :: Rj(3), Rk(3), Fi(3), Fk(3), Fl(3), Fij(3)
    real(rb)    :: normRkj, normX, a, b, phi
    real(rb)    :: rij(3), rkj(3), rlk(3), x(3), y(3), z(3), u(3), v(3), w(3)

    integer(ib),     pointer :: type(:)
    real(rb),        pointer :: R(:,:), F(:,:), charge(:)
    type(tModel),    pointer :: model
    type(tStruct),   pointer :: dihedral(:), d
    type(tModelPtr), pointer :: pairType(:,:)

    call c_f_pointer( me%dihedral%list, dihedral, [me%dihedral%number] )
    call c_f_pointer( me%R, R, [3,me%natoms] )
    call c_f_pointer( me%F, F, [3,me%natoms] )
    call c_f_pointer( me%charge, charge, [me%natoms] )
    call c_f_pointer( me%type, type, [me%natoms] )
    call c_f_pointer( me%pairType, pairType, [me%ntypes,me%ntypes] )

    invL = one/L
    Energy = zero
    Virial = zero
    do m = 1, me%dihedral%number
      d => dihedral(m)
      Rj = R(:,d%j)
      Rk = R(:,d%k)
      rij = R(:,d%i) - Rj
      rkj = Rk - Rj
      rlk = R(:,d%l) - Rk
      rij = rij - L*nint(rij*InvL)
      rkj = rkj - L*nint(rkj*InvL)
      rlk = rlk - L*nint(rlk*InvL)
      normRkj = sqrt(sum(rkj*rkj))
      z = rkj/normRkj
      x = rij - sum(rij*z)*z
      normX = sqrt(sum(x*x))
      x = x/normX
      y = cross(z,x)
      a = sum(x*rlk)
      b = sum(y*rlk)
      phi = atan2(b,a)
      model => d%model
      call compute_dihedral()
      Fd = Fd/(a*a + b*b)
      u = (a*cross(rlk,z) - b*rlk)/normX
      v = (a*cross(rlk,x) + sum(z*u)*rij)/normRkj
      w = v + sum(z*rij)*u/normRkj
      Fi = Fd*sum(u*y)*y
      Fl = Fd*(a*y - b*x)
      Fk = -(Fd*(sum(v*x)*x + sum(w*y)*y) + Fl)
      F(:,d%i) = F(:,d%i) + Fi
      F(:,d%k) = F(:,d%k) + Fk
      F(:,d%l) = F(:,d%l) + Fl
      F(:,d%j) = F(:,d%j) + (Fi + Fk + Fl)
      Energy = Energy + Ed
      Virial = Virial + sum(Fi*rij + Fk*rkj + Fl*(rlk + rkj))
      if (model%f14 /= zero) then
        rij = rij + rlk - rkj
        r2 = sum(rij*rij)
        if (r2 < me%RcSq) then
          invR2 = one/r2
          model => pairType(type(d%i),type(d%l))%model
          icharge = charge(d%i)
          jcharge = charge(d%l)
          call compute_pair()
          Eij = model%f14*Eij
          Wij = model%f14*Wij
          Energy = Energy + Eij
          Virial = Virial + Wij
          Fij = Wij*invR2*rij
          F(:,d%i) = F(:,d%i) + Fij
          F(:,d%l) = F(:,d%l) - Fij
        end if
      end if
    end do
    me%dihedral%energy = Energy
    me%dihedral%virial = third*Virial
    nullify( R, F, charge, type, pairType, dihedral )
    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      function cross( a, b ) result( c )
        real(rb), intent(in) :: a(3), b(3)
        real(rb) :: c(3)
        c = [ a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1) ]
      end function cross
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      include "compute_dihedral.f90"
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      include "compute_pair.f90"
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine compute_dihedrals

!---------------------------------------------------------------------------------------------------

  subroutine find_pairs_and_compute( me, L )
    type(tEmDee), intent(inout) :: me
    real(rb),     intent(in)    :: L

    integer(ib) :: i, j, k, m, n, icell, jcell, maxpairs, npairs, nlocal, kprev, mprev, itype
    real(rb)    :: invL, invL2, xRc2, Rc2, r2, invR2, Epot, Virial, Eij, Wij, icharge, jcharge
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
                jcharge = charge(j)
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
    real(rb) :: invL, invL2, Rc2, r2, invR2, Epot, Virial, Eij, Wij, icharge, jcharge
    real(rb) :: Rs(3,me%natoms), Rij(3), Ri(3), Fi(3), Fij(3)

    type(tModel),    pointer :: model
    integer(ib),     pointer :: type(:), first(:), last(:), neighbor(:)
    real(rb),        pointer :: R(:,:), F(:,:), charge(:)
    type(tModelPtr), pointer :: pairType(:,:)

    real(rb), allocatable :: Fj(:,:)

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
          jcharge = charge(j)
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
