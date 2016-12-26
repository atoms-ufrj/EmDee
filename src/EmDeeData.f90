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

module EmDeeData

use omp_lib
use lists
use math
use models
use structs
use ArBee

implicit none

integer, parameter :: extra = 2000

integer, parameter :: ndiv = 2
integer, parameter :: nbcells = 62
integer, parameter :: nb(3,nbcells) = reshape( [ &
   1, 0, 0,    2, 0, 0,   -2, 1, 0,   -1, 1, 0,    0, 1, 0,    1, 1, 0,    2, 1, 0,   -2, 2, 0,  &
  -1, 2, 0,    0, 2, 0,    1, 2, 0,    2, 2, 0,   -2,-2, 1,   -1,-2, 1,    0,-2, 1,    1,-2, 1,  &
   2,-2, 1,   -2,-1, 1,   -1,-1, 1,    0,-1, 1,    1,-1, 1,    2,-1, 1,   -2, 0, 1,   -1, 0, 1,  &
   0, 0, 1,    1, 0, 1,    2, 0, 1,   -2, 1, 1,   -1, 1, 1,    0, 1, 1,    1, 1, 1,    2, 1, 1,  &
  -2, 2, 1,   -1, 2, 1,    0, 2, 1,    1, 2, 1,    2, 2, 1,   -2,-2, 2,   -1,-2, 2,    0,-2, 2,  &
   1,-2, 2,    2,-2, 2,   -2,-1, 2,   -1,-1, 2,    0,-1, 2,    1,-1, 2,    2,-1, 2,   -2, 0, 2,  &
  -1, 0, 2,    0, 0, 2,    1, 0, 2,    2, 0, 2,   -2, 1, 2,   -1, 1, 2,    0, 1, 2,    1, 1, 2,  &
   2, 1, 2,   -2, 2, 2,   -1, 2, 2,    0, 2, 2,    1, 2, 2,    2, 2, 2 ], shape(nb) )

type, private :: tCell
  integer :: neighbor(nbcells)
end type tCell

type :: tData

  integer :: natoms                       ! Number of atoms in the system
  integer :: mcells = 0                   ! Number of cells at each dimension
  integer :: ncells = 0                   ! Total number of cells
  integer :: maxcells = 0                 ! Maximum number of cells
  integer :: maxatoms = 0                 ! Maximum number of atoms in a cell
  integer :: maxpairs = 0                 ! Maximum number of pairs formed by all atoms of a cell
  integer :: ntypes                       ! Number of atom types
  integer :: nbodies = 0                  ! Number of rigid bodies
  integer :: maxbodies = 0                ! Maximum number of rigid bodies
  integer :: nfree                        ! Number of independent atoms
  integer :: nthreads                     ! Number of parallel openmp threads
  integer :: threadAtoms                  ! Number of atoms per parallel thread
  integer :: threadBodies = 0             ! Number of rigid bodies per parallel thread
  integer :: nlayers = 1                  ! Number of pair interaction model layers
  integer :: layer = 1                    ! Current layer of pair interaction models

  real(rb) :: Lbox = zero                 ! Length of the simulation box
  real(rb) :: invL = huge(one)            ! Inverse length of the simulation box
  real(rb) :: invL2 = huge(one)           ! Squared inverse length of the simulation box

  real(rb) :: Rc                          ! Cut-off distance
  real(rb) :: RcSq                        ! Cut-off distance squared
  real(rb) :: xRc                         ! Extended cutoff distance (including skin)
  real(rb) :: xRcSq                       ! Extended cutoff distance squared
  real(rb) :: skinSq                      ! Square of the neighbor list skin width
  real(rb) :: totalMass                   ! Sum of the masses of all atoms
  real(rb) :: startTime                   ! Time recorded at initialization
  real(rb) :: eshift                      ! Potential shifting factor for Coulombic interactions
  real(rb) :: fshift                      ! Force shifting factor for Coulombic interactions

  logical :: initialized = .false.        ! Flag for checking coordinates initialization

  type(kiss) :: random                    ! Random number generator

  type(tList) :: cellAtom                 ! List of atoms belonging to each cell
  type(tList) :: threadCell               ! List of cells to be dealt with in each parallel thread
  type(tList) :: excluded                 ! List of pairs excluded from the neighbor lists

  type(structList) :: bonds               ! List of bonds
  type(structList) :: angles              ! List of angles
  type(structList) :: dihedrals           ! List of dihedrals

  integer, allocatable :: atomType(:)     ! Type indexes of all atoms
  integer, allocatable :: atomCell(:)     ! Array containing the current cell of each atom
  integer, allocatable :: free(:)         ! Pointer to the list of independent atoms

  real(rb), allocatable :: R(:,:)         ! Coordinates of all atoms
  real(rb), allocatable :: P(:,:)         ! Momenta of all atoms
  real(rb), allocatable :: F(:,:)         ! Resultant forces on all atoms
  real(rb), allocatable :: charge(:)      ! Electric charges of all atoms
  real(rb), allocatable :: mass(:)        ! Masses of all atoms
  real(rb), allocatable :: invMass(:)     ! Inverses of atoms masses
  real(rb), allocatable :: R0(:,:)        ! Position of each atom at latest neighbor list building

  type(tCell), allocatable :: cell(:)      ! Array containing all neighbor cells of each cell
  type(tBody), allocatable :: body(:)      ! Pointer to the rigid bodies present in the system
  type(tList), allocatable :: neighbor(:)  ! Pointer to neighbor lists

  type(pairModelContainer), allocatable :: pair(:,:,:)
  logical,  allocatable :: multilayer(:,:)
  logical,  allocatable :: overridable(:,:)
  real(rb), pointer     :: layer_energy(:)

end type tData

contains

!===================================================================================================

  subroutine set_pair_type( me, itype, jtype, layer, container )
    type(tData),          intent(inout) :: me
    integer(ib),          intent(in)    :: itype, jtype, layer
    type(modelContainer), intent(in)    :: container

    integer :: ktype

    select type (pmodel => container%model)
      class is (cPairModel)
        associate (pair => me%pair(:,:,layer))
          if (itype == jtype) then
            pair(itype,itype) = container
            call pair(itype,itype) % model % shifting_setup( me%Rc )
            do ktype = 1, me%ntypes
              if ((ktype /= itype).and.me%overridable(itype,ktype)) then
                pair(itype,ktype) = pair(ktype,ktype) % mix( pmodel )
                call pair(itype,ktype) % model % shifting_setup( me%Rc )
                pair(ktype,itype) = pair(itype,ktype)
              end if
            end do
          else
            pair(itype,jtype) = container
            call pair(itype,jtype) % model % shifting_setup( me%Rc )
            pair(jtype,itype) = pair(itype,jtype)
          end if
        end associate
      class default
        stop "ERROR: a valid pair model must be provided"
    end select

  end subroutine set_pair_type

!===================================================================================================

  subroutine rigid_body_forces( me, Virial )
    type(tData), intent(inout) :: me
    real(rb),    intent(inout) :: Virial

    real(rb) :: Wrb(me%nthreads)

    !$omp parallel num_threads(me%nthreads)
    block
      integer :: thread
      thread = omp_get_thread_num() + 1
      call compute_body_forces( thread, Wrb(thread) )
    end block
    !$omp end parallel
    Virial = Virial - third*sum(Wrb)

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine compute_body_forces( thread, Wrb )
        integer,  intent(in)  :: thread
        real(rb), intent(out) :: Wrb
        integer :: i
        Wrb = zero
        do i = (thread - 1)*me%threadBodies + 1, min(thread*me%threadBodies, me%nbodies)
          associate (b => me%body(i))
            call b % force_torque_virial( me%F )
            Wrb = Wrb + b%virial
          end associate
        end do
      end subroutine compute_body_forces
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine rigid_body_forces

!===================================================================================================

  real(rb) function maximum_approach_sq( N, delta )
    integer,  intent(in) :: N
    real(rb), intent(in) :: delta(3,N)

    integer  :: i
    real(rb) :: maximum, next, deltaSq

    maximum = sum(delta(:,1)**2)
    next = maximum
    do i = 2, N
      deltaSq = sum(delta(:,i)**2)
      if (deltaSq > maximum) then
        next = maximum
        maximum = deltaSq
      end if
    end do
    maximum_approach_sq = maximum + 2*sqrt(maximum*next) + next

  end function maximum_approach_sq

!===================================================================================================

  subroutine distribute_atoms( me, M, Rs )
    type(tData), intent(inout) :: me
    integer, intent(in)    :: M
    real(rb),    intent(in)    :: Rs(3,me%natoms)

    integer :: MM, cells_per_thread, maxNatoms, threadNatoms(me%nthreads), next(me%natoms)
    logical :: make_cells
    integer, allocatable :: natoms(:)

    MM = M*M
    make_cells = M /= me%mcells
    if (make_cells) then
      me%mcells = M
      me%ncells = M*MM
      if (me%ncells > me%maxcells) then
        deallocate( me%cell, me%cellAtom%first, me%cellAtom%last )
        allocate( me%cell(me%ncells), me%cellAtom%first(me%ncells), me%cellAtom%last(me%ncells) )
        call me % threadCell % allocate( 0, me%nthreads )
        me%maxcells = me%ncells
      end if
      cells_per_thread = (me%ncells + me%nthreads - 1)/me%nthreads
    end if

    allocate( natoms(me%ncells) )

    !$omp parallel num_threads(me%nthreads) reduction(max:maxNatoms)
    call distribute( omp_get_thread_num() + 1, maxNatoms )
    !$omp end parallel
    me%maxatoms = maxNatoms
    me%maxpairs = (maxNatoms*((2*nbcells + 1)*maxNatoms - 1))/2

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine distribute( thread, maxNatoms )
        integer, intent(in)  :: thread
        integer, intent(out) :: maxNatoms

        integer :: i, k, icell, ix, iy, iz, first, last, atoms_per_thread
        integer :: icoord(3)
        integer, allocatable :: head(:)

        if (make_cells) then
          first = (thread - 1)*cells_per_thread + 1
          last = min( thread*cells_per_thread, me%ncells )
          do icell = first, last
            k = icell - 1
            iz = k/MM
            iy = (k - iz*MM)/M
            ix = k - (iy*M + iz*MM)
            me%cell(icell)%neighbor = 1 + pbc(ix+nb(1,:)) + pbc(iy+nb(2,:))*M + pbc(iz+nb(3,:))*MM
          end do
          me%threadCell%first(thread) = first
          me%threadCell%last(thread) = last
        else
          first = me%threadCell%first(thread)
          last = me%threadCell%last(thread)
        end if

        atoms_per_thread = (me%natoms + me%nthreads - 1)/me%nthreads
        do i = (thread - 1)*atoms_per_thread + 1, min( thread*atoms_per_thread, me%natoms )
          icoord = int(M*(Rs(:,i) - floor(Rs(:,i))),ib)
          me%atomCell(i) = 1 + icoord(1) + M*icoord(2) + MM*icoord(3)
        end do
        !$omp barrier

        allocate( head(first:last) )
        head = 0
        natoms(first:last) = 0
        do i = 1, me%natoms
          icell = me%atomCell(i)
          if ((icell >= first).and.(icell <= last)) then
            next(i) = head(icell)
            head(icell) = i
            natoms(icell) = natoms(icell) + 1
          end if
        end do
        threadNatoms(thread) = sum(natoms(first:last))
        !$omp barrier

        maxNatoms = 0
        k = sum(threadNatoms(1:thread-1))
        do icell = first, last
          me%cellAtom%first(icell) = k + 1
          i = head(icell)
          do while (i /= 0)
            k = k + 1
            me%cellAtom%item(k) = i
            i = next(i)
          end do
          me%cellAtom%last(icell) = k
          if (natoms(icell) > maxNatoms) maxNatoms = natoms(icell)
        end do
      end subroutine distribute
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      elemental integer function pbc( x )
        integer, intent(in) :: x
        if (x < 0) then
          pbc = x + M
        else if (x >= M) then
          pbc = x - M
        else
          pbc = x
        end if
      end function pbc
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine distribute_atoms

!===================================================================================================

  subroutine compute_bonds( me, threadId, R, F, Potential, Virial )
    type(tData), intent(inout) :: me
    integer,     intent(in)    :: threadId
    real(rb),    intent(in)    :: R(3,me%natoms)
    real(rb),    intent(inout) :: F(3,me%natoms), Potential, Virial

    integer  :: m, nbonds
    real(rb) :: invR2, E, W
    real(rb) :: Rij(3), Fij(3)

    nbonds = (me%bonds%number + me%nthreads - 1)/me%nthreads
    do m = (threadId - 1)*nbonds + 1, min( me%bonds%number, threadId*nbonds )
      associate(bond => me%bonds%item(m))
        Rij = R(:,bond%i) - R(:,bond%j)
        Rij = Rij - anint(Rij)
        invR2 = me%invL2/sum(Rij*Rij)
        select type (model => bond%model)
          include "compute_bond.f90"
        end select
        Potential = Potential + E
        Virial = Virial + W
        Fij = W*Rij*invR2
        F(:,bond%i) = F(:,bond%i) + Fij
        F(:,bond%j) = F(:,bond%j) - Fij
      end associate
    end do

  end subroutine compute_bonds

!===================================================================================================

  subroutine compute_angles( me, threadId, R, F, Potential, Virial )
    type(tData), intent(inout) :: me
    integer,     intent(in)    :: threadId
    real(rb),    intent(in)    :: R(3,me%natoms)
    real(rb),    intent(inout) :: F(3,me%natoms), Potential, Virial

    integer  :: i, j, k, m, nangles
    real(rb) :: aa, bb, ab, axb, theta, Ea, Fa
    real(rb) :: Rj(3), Fi(3), Fk(3), a(3), b(3)

    nangles = (me%angles%number + me%nthreads - 1)/me%nthreads
    do m = (threadId - 1)*nangles + 1, min( me%angles%number, threadId*nangles )
      associate (angle => me%angles%item(m))
        i = angle%i
        j = angle%j
        k = angle%k
        Rj = R(:,j)
        a = R(:,i) - Rj
        b = R(:,k) - Rj
        a = a - anint(a)
        b = b - anint(b)
        aa = sum(a*a)
        bb = sum(b*b)
        ab = sum(a*b)
        axb = sqrt(aa*bb - ab*ab)
        theta = atan2(axb,ab)
        select type (model => angle%model)
          include "compute_angle.f90"
        end select
        Fa = Fa/(me%Lbox*axb)
        Fi = Fa*(b - (ab/aa)*a)
        Fk = Fa*(a - (ab/bb)*b)
        F(:,i) = F(:,i) + Fi
        F(:,k) = F(:,k) + Fk
        F(:,j) = F(:,j) - (Fi + Fk)
        Potential = Potential + Ea
        Virial = Virial + me%Lbox*sum(Fi*a + Fk*b)
      end associate
    end do

  end subroutine compute_angles

!===================================================================================================

  subroutine compute_dihedrals( me, threadId, R, F, Potential, Virial )
    type(tData), intent(inout) :: me
    integer,     intent(in)    :: threadId
    real(rb),    intent(in)    :: R(3,me%natoms)
    real(rb),    intent(inout) :: F(3,me%natoms), Potential, Virial

    integer  :: m, ndihedrals, i, j
    real(rb) :: Rc2, Ed, Fd, r2, invR2, Eij, Wij, Qi, Qj
    real(rb) :: Rj(3), Rk(3), Fi(3), Fk(3), Fl(3), Fij(3)
    real(rb) :: normRkj, normX, a, b, phi, factor14
    real(rb) :: rij(3), rkj(3), rlk(3), x(3), y(3), z(3), u(3), v(3), w(3)

    Rc2 = me%RcSq*me%invL2
    ndihedrals = (me%dihedrals%number + me%nthreads - 1)/me%nthreads
    do m = (threadId - 1)*ndihedrals + 1, min( me%dihedrals%number, threadId*ndihedrals )
      associate (dihedral => me%dihedrals%item(m))
        Rj = R(:,dihedral%j)
        Rk = R(:,dihedral%k)
        rij = R(:,dihedral%i) - Rj
        rkj = Rk - Rj
        rlk = R(:,dihedral%l) - Rk
        rij = rij - anint(rij)
        rkj = rkj - anint(rkj)
        rlk = rlk - anint(rlk)
        normRkj = sqrt(sum(rkj*rkj))
        z = rkj/normRkj
        x = rij - sum(rij*z)*z
        normX = sqrt(sum(x*x))
        x = x/normX
        y = cross(z,x)
        a = sum(x*rlk)
        b = sum(y*rlk)
        phi = atan2(b,a)
        select type (model => dihedral%model)
          class is (cDihedralModel)
            factor14 = model%factor14
          include "compute_dihedral.f90"
        end select
        Fd = Fd/(me%Lbox*(a*a + b*b))
        u = (a*cross(rlk,z) - b*rlk)/normX
        v = (a*cross(rlk,x) + sum(z*u)*rij)/normRkj
        w = v + sum(z*rij)*u/normRkj
        Fi = Fd*sum(u*y)*y
        Fl = Fd*(a*y - b*x)
        Fk = -(Fd*(sum(v*x)*x + sum(w*y)*y) + Fl)
        F(:,dihedral%i) = F(:,dihedral%i) + Fi
        F(:,dihedral%k) = F(:,dihedral%k) + Fk
        F(:,dihedral%l) = F(:,dihedral%l) + Fl
        F(:,dihedral%j) = F(:,dihedral%j) + (Fi + Fk + Fl)
        Potential = Potential + Ed
        Virial = Virial + me%Lbox*sum(Fi*rij + Fk*rkj + Fl*(rlk + rkj))
        if (factor14 /= zero) then
          i = dihedral%i
          j = dihedral%l
          rij = rij + rlk - rkj
          r2 = sum(rij*rij)
          if (r2 < me%RcSq) then
            invR2 = me%invL2/r2
            Qi = me%charge(i)
            Qj = me%charge(j)
            select type ( model => me%pair(me%atomType(i),me%atomType(j),me%layer)%model )
              include "compute_pair.f90"
            end select
            Eij = factor14*Eij
            Wij = factor14*Wij
            Potential = Potential + Eij
            Virial = Virial + Wij
            Fij = Wij*invR2*rij
            F(:,i) = F(:,i) + Fij
            F(:,j) = F(:,j) - Fij
          end if
        end if
      end associate
    end do

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      function cross( a, b ) result( c )
        real(rb), intent(in) :: a(3), b(3)
        real(rb) :: c(3)
        c = [ a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1) ]
      end function cross
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine compute_dihedrals

!===================================================================================================

  subroutine find_pairs_and_compute( me, thread, Rs, F, Potential, Virial, Elayer )
    type(tData), intent(inout) :: me
    integer,     intent(in)    :: thread
    real(rb),    intent(in)    :: Rs(3,me%natoms)
    real(rb),    intent(out)   :: F(3,me%natoms), Potential, Virial, Elayer(me%nlayers)

    integer  :: i, j, k, m, n, icell, jcell, npairs, itype, jtype, layer
    integer  :: nlocal, ntotal, first, last
    real(rb) :: xRc2, Rc2, r2, invR2, Eij, Wij, Qi, Qj
    logical  :: include(0:me%maxpairs)
    integer  :: atom(me%maxpairs), index(me%natoms)
    real(rb) :: Ri(3), Rij(3), Fi(3), Fij(3)
    integer,  allocatable :: xlist(:)

    xRc2 = me%xRcSq*me%invL2
    Rc2 = me%RcSq*me%invL2

    F = zero
    Potential = zero
    Virial = zero
    Elayer = zero

    include = .true.
    index = 0
    npairs = 0
    associate (neighbor => me%neighbor(thread))
      do icell = me%threadCell%first(thread), me%threadCell%last(thread)

        if (neighbor%nitems < npairs + me%maxpairs) then
          call neighbor % resize( npairs + me%maxpairs + extra )
        end if

        first = me%cellAtom%first(icell)
        last = me%cellAtom%last(icell)
        nlocal = last - first + 1
        atom(1:nlocal) = me%cellAtom%item(first:last)

        ntotal = nlocal
        do m = 1, nbcells
          jcell = me%cell(icell)%neighbor(m)
          first = me%cellAtom%first(jcell)
          last = me%cellAtom%last(jcell)
          n = ntotal + 1
          ntotal = n + last - first
          atom(n:ntotal) = me%cellAtom%item(first:last)
        end do

        forall (m=1:ntotal) index(atom(m)) = m
        do k = 1, nlocal
          i = atom(k)
          neighbor%first(i) = npairs + 1
          itype = me%atomType(i)
          Qi = me%charge(i)
          Ri = Rs(:,i)
          Fi = zero
          xlist = index(me%excluded%item(me%excluded%first(i):me%excluded%last(i)))
          include(xlist) = .false.
          associate (partner => me%pair(:,itype,me%layer), multilayer => me%multilayer(:,itype))
            do m = k + 1, ntotal
              if (include(m)) then
                j = atom(m)
                jtype = me%atomType(j)
                select type ( model => partner(jtype)%model )
                  type is (pair_none)
                  class default
                    Qj = me%charge(j)
                    Rij = Ri - Rs(:,j)
                    Rij = Rij - anint(Rij)
                    r2 = sum(Rij*Rij)
                    if (r2 < xRc2) then
                      npairs = npairs + 1
                      neighbor%item(npairs) = j
                      if (r2 < Rc2) then
                        invR2 = me%invL2/r2
                        select type (model)
                          include "compute_pair.f90"
                        end select
                        Potential = Potential + Eij
                        Virial = Virial + Wij
                        Fij = Wij*invR2*Rij
                        Fi = Fi + Fij
                        F(:,j) = F(:,j) - Fij

                        if (multilayer(jtype)) then
                          do layer = 1, me%nlayers
                            select type ( model => me%pair(itype,jtype,layer)%model )
                              include "compute_pair.f90"
                            end select
                            Elayer(layer) = Elayer(layer) + Eij
                          end do
                        end if

                      end if
                    end if
                end select
              end if
            end do
          end associate
          F(:,i) = F(:,i) + Fi
          neighbor%last(i) = npairs
          include(xlist) = .true.
        end do
        index(atom(1:ntotal)) = 0

      end do
      neighbor%count = npairs
    end associate

  end subroutine find_pairs_and_compute

!===================================================================================================

  subroutine compute_pairs( me, thread, Rs, F, Potential, Virial, Elayer )
    type(tData), intent(in)  :: me
    integer,     intent(in)  :: thread
    real(rb),    intent(in)  :: Rs(3,me%natoms)
    real(rb),    intent(out) :: F(3,me%natoms), Potential, Virial, Elayer(me%nlayers)

    integer  :: i, j, k, m, itype, jtype, firstAtom, lastAtom, layer
    real(rb) :: Rc2, r2, invR2, Eij, Wij, Qi, Qj
    real(rb) :: Rij(3), Ri(3), Fi(3), Fij(3)

    Rc2 = me%RcSq*me%invL2

    F = zero
    Potential = zero
    Virial = zero
    Elayer = zero

    firstAtom = me%cellAtom%first(me%threadCell%first(thread))
    lastAtom = me%cellAtom%last(me%threadCell%last(thread))
    associate (neighbor => me%neighbor(thread), Q => me%charge)
      do m = firstAtom, lastAtom
        i = me%cellAtom%item(m)
        itype = me%atomType(i)
        Qi = Q(i)
        Ri = Rs(:,i)
        Fi = zero
        associate (pair => me%pair(:,itype,me%layer), multilayer => me%multilayer(:,itype))
          do k = neighbor%first(i), neighbor%last(i)
            j = neighbor%item(k)
            Qj = Q(j)
            Rij = Ri - Rs(:,j)
            Rij = Rij - anint(Rij)
            r2 = sum(Rij*Rij)
            if (r2 < Rc2) then
              invR2 = me%invL2/r2
              jtype = me%atomType(j)
              select type ( model => pair(jtype)%model )
                include "compute_pair.f90"
              end select
              Potential = Potential + Eij
              Virial = Virial + Wij
              Fij = Wij*invR2*Rij
              Fi = Fi + Fij
              F(:,j) = F(:,j) - Fij

              if (multilayer(jtype)) then
                do layer = 1, me%nlayers
                  select type ( model => me%pair(itype,jtype,layer)%model )
                    include "compute_pair.f90"
                  end select
                  Elayer(layer) = Elayer(layer) + Eij
                end do
              end if

            end if
          end do
        end associate
        F(:,i) = F(:,i) + Fi
      end do
    end associate

  end subroutine compute_pairs

!===================================================================================================

  subroutine update_pairs( me, thread, Rs, DF, DE, DW, new_layer )
    type(tData), intent(in)  :: me
    integer,     intent(in)  :: thread, new_layer
    real(rb),    intent(in)  :: Rs(3,me%natoms)
    real(rb),    intent(out) :: DF(3,me%natoms), DE, DW

    integer  :: i, j, k, m, itype, jtype, firstAtom, lastAtom
    real(rb) :: Rc2, r2, invR2, new_Eij, new_Wij, Eij, Wij, Qi, Qj
    real(rb) :: Rij(3), Ri(3), Fi(3), Fij(3)
    logical  :: participant(me%ntypes)

    participant = any(me%multilayer,dim=2)
    Rc2 = me%RcSq*me%invL2

    DF = zero
    DE = zero
    DW = zero

    firstAtom = me%cellAtom%first(me%threadCell%first(thread))
    lastAtom = me%cellAtom%last(me%threadCell%last(thread))
    associate (neighbor => me%neighbor(thread), Q => me%charge)
      do m = firstAtom, lastAtom
        i = me%cellAtom%item(m)
        itype = me%atomType(i)
        if (participant(itype)) then
          Qi = Q(i)
          Ri = Rs(:,i)
          Fi = zero
          associate (multilayer => me%multilayer(:,itype))
            do k = neighbor%first(i), neighbor%last(i)
              j = neighbor%item(k)
              jtype = me%atomType(j)
              if (multilayer(jtype)) then
                Qj = Q(j)
                Rij = Ri - Rs(:,j)
                Rij = Rij - anint(Rij)
                r2 = sum(Rij*Rij)
                if (r2 < Rc2) then
                  invR2 = me%invL2/r2
                  select type ( model => me%pair(jtype,itype,new_layer)%model )
                    include "compute_pair.f90"
                  end select
                  new_Eij = Eij
                  new_Wij = Wij
                  select type ( model => me%pair(jtype,itype,me%layer)%model )
                    include "compute_pair.f90"
                  end select
                  DE = DE + new_Eij - Eij
                  Wij = new_Wij - Wij
                  DW = DW + Wij
                  Fij = Wij*invR2*Rij
                  Fi = Fi + Fij
                  DF(:,j) = DF(:,j) - Fij
                end if
              end if
            end do
          end associate
          DF(:,i) = DF(:,i) + Fi
        end if
      end do
    end associate

  end subroutine update_pairs

!===================================================================================================

end module EmDeeData
