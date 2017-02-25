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
 0,  0,  1,    0,  1,  0,    1,  0,  0,   -1,  0,  1,   -1,  1,  0,    0, -1,  1,    0,  1,  1, &
 1,  0,  1,    1,  1,  0,   -1, -1,  1,   -1,  1,  1,    1, -1,  1,    1,  1,  1,    0,  0,  2, &
 0,  2,  0,    2,  0,  0,   -2,  0,  1,   -2,  1,  0,   -1,  0,  2,   -1,  2,  0,    0, -2,  1, &
 0, -1,  2,    0,  1,  2,    0,  2,  1,    1,  0,  2,    1,  2,  0,    2,  0,  1,    2,  1,  0, &
-2, -1,  1,   -2,  1,  1,   -1, -2,  1,   -1, -1,  2,   -1,  1,  2,   -1,  2,  1,    1, -2,  1, &
 1, -1,  2,    1,  1,  2,    1,  2,  1,    2, -1,  1,    2,  1,  1,   -2,  0,  2,   -2,  2,  0, &
 0, -2,  2,    0,  2,  2,    2,  0,  2,    2,  2,  0,   -2, -2,  1,   -2, -1,  2,   -2,  1,  2, &
-2,  2,  1,   -1, -2,  2,   -1,  2,  2,    1, -2,  2,    1,  2,  2,    2, -2,  1,    2, -1,  2, &
 2,  1,  2,    2,  2,  1,   -2, -2,  2,   -2,  2,  2,    2, -2,  2,    2,  2,  2], shape(nb) )

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
  integer :: threadFreeAtoms              ! Number of free atoms per parallel thread
  integer :: threadBodies = 0             ! Number of rigid bodies per parallel thread
  integer :: nlayers = 1                  ! Number of pair interaction model layers
  integer :: layer = 1                    ! Current layer of pair interaction models

  real(rb) :: Lbox = zero                 ! Length of the simulation box
  real(rb) :: Lbox3(3) = zero
  real(rb) :: invL = huge(one)            ! Inverse length of the simulation box
  real(rb) :: invL2 = huge(one)           ! Squared inverse length of the simulation box

  real(rb) :: Rc                          ! Cut-off distance
  real(rb) :: skin                        ! Neighbor list skin width
  real(rb) :: RcSq                        ! Cut-off distance squared
  real(rb) :: xRc                         ! Extended cutoff distance (including skin)
  real(rb) :: xRcSq                       ! Extended cutoff distance squared
  real(rb) :: skinSq                      ! Square of the neighbor list skin width
  real(rb) :: totalMass                   ! Sum of the masses of all atoms
  real(rb) :: startTime                   ! Time recorded at initialization

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
  integer, allocatable :: atomBody(:)     ! Array containing the rigid body containing each atom
  integer, allocatable :: free(:)         ! Pointer to the list of independent atoms

  real(rb), allocatable :: R(:,:)         ! Coordinates of all atoms
  real(rb), allocatable :: P(:,:)         ! Momenta of all atoms
  real(rb), allocatable :: F(:,:)         ! Resultant forces on all atoms
  real(rb), allocatable :: charge(:)      ! Electric charges of all atoms
  real(rb), allocatable :: mass(:)        ! Masses of all atoms
  real(rb), allocatable :: invMass(:)     ! Inverses of atoms masses
  real(rb), allocatable :: R0(:,:)        ! Position of each atom at latest neighbor list building
  logical,  allocatable :: charged(:)     ! Flag to determine if a particle is charged

  type(tCell), allocatable :: cell(:)      ! Array containing all neighbor cells of each cell
  type(tBody), allocatable :: body(:)      ! Pointer to the rigid bodies present in the system
  type(tList), allocatable :: neighbor(:)  ! Pointer to neighbor lists

  type(pairModelContainer), allocatable :: pair(:,:,:)
  type(coulModelContainer), allocatable :: coul(:)
  class(cKspaceModel),      allocatable :: kspace

  logical,  allocatable :: multilayer(:,:)
  logical,  allocatable :: overridable(:,:)
  logical,  allocatable :: interact(:,:)
  integer,  allocatable :: other_layer(:)
  real(rb), allocatable :: threadEnergy(:,:)
  real(rb), pointer     :: Energy(:)

  logical :: multilayer_coulomb
  logical :: kspace_active

  logical     :: respa                    ! Flag to determine if RESPA was activated
  real(rb)    :: respaRc                   ! Internal cutoff distance for RESPA
  real(rb)    :: respaRcSq                 ! Internal cutoff distance squared
  real(rb)    :: xRespaRcSq                ! Extended internal cutoff distance squared
  real(rb)    :: fshift
  integer(ib) :: Npair                    ! Number of core-neighbor sweepings per MD step
  integer(ib) :: Nbond                    ! Number of bonded computations per core sweeping

  type(pairModelContainer), allocatable :: shortPair(:,:,:)
  real(rb),                 allocatable :: shortF(:,:)

end type tData

contains

!===================================================================================================

  subroutine assign_momenta( me, thread, P, twoKEt, twoKEr )
    type(tData), intent(inout) :: me
    integer,     intent(in)    :: thread
    real(rb),    intent(in)    :: P(3,me%natoms)
    real(rb),    intent(out)   :: twoKEt(3), twoKEr(3)

    integer :: i, j
    real(rb) :: L(3), Pj(3)

    twoKEt = zero
    twoKEr = zero
    do j = (thread - 1)*me%threadFreeAtoms + 1, min(thread*me%threadFreeAtoms, me%nfree)
      i = me%free(j)
      me%P(:,i) = P(:,i)
      twoKEt = twoKEt + me%invMass(i)*P(:,i)**2
    end do
    do i = (thread - 1)*me%threadBodies + 1, min(thread*me%threadBodies, me%nbodies)
      associate(b => me%body(i))
        b%pcm = zero
        L = zero
        do j = 1, b%NP
          Pj = P(:,b%index(j))
          b%pcm = b%pcm + Pj
          L = L + cross_product( b%delta(:,j), Pj )
        end do
        twoKEt = twoKEt + b%invMass*b%pcm**2
        b%pi = matmul( matrix_C(b%q), two*L )
        b%omega= half*b%invMoI*matmul( matrix_Bt(b%q), b%pi )
        twoKEr = twoKEr + b%MoI*b%omega**2
      end associate
    end do

  end subroutine assign_momenta

!===================================================================================================

  subroutine check_actual_interactions( me )
    type(tData), intent(inout) :: me

    integer :: i, j, k
    logical :: no_pair, no_coul, coul_only, neutral(me%ntypes), inert
    type(pair_none) :: PairNone
    type(coul_none) :: CoulNone

    do i = 1, me%ntypes
      neutral(i) = count((me%atomType == i) .and. me%charged) == 0
      do j = 1, i
        associate( pair => me%pair(i,j,:) )
          k = 0
          inert = .true.
          do while (inert .and. (k < me%nlayers))
            k = k + 1
            no_pair = same_type_as( pair(k)%model, PairNone )
            no_coul = same_type_as( me%coul(k)%model, CoulNone ) .or. (.not.pair(k)%coulomb)
            coul_only = no_pair .and. (.not.no_coul)
            inert = (no_pair .and. no_coul) .or. (coul_only .and. neutral(i) .and. neutral(j))
          end do
          me%interact(i,j) = .not.inert
          me%interact(j,i) = me%interact(i,j)
        end associate
      end do
    end do

  end subroutine check_actual_interactions

!===================================================================================================

  subroutine set_pair_type( me, itype, jtype, layer, container, kCoul )
    type(tData),          intent(inout) :: me
    integer(ib),          intent(in)    :: itype, jtype, layer
    type(modelContainer), intent(in)    :: container
    real(rb),             intent(in)    :: kCoul

    integer :: ktype

    select type (pmodel => container%model)
      class is (cPairModel)
        associate (pair => me%pair(:,:,layer))
          if (itype == jtype) then
            pair(itype,itype) = container
            pair(itype,itype)%coulomb = kCoul /= zero
            if (pair(itype,itype)%coulomb) pair(itype,itype)%kCoul = kCoul
            call pair(itype,itype) % model % shifting_setup( me%Rc )
            do ktype = 1, me%ntypes
              if ((ktype /= itype).and.me%overridable(itype,ktype)) then
                pair(itype,ktype) = pair(ktype,ktype) % mix( pair(itype,itype) )
                call pair(itype,ktype) % model % shifting_setup( me%Rc )
                pair(ktype,itype) = pair(itype,ktype)
              end if
            end do
          else
            pair(itype,jtype) = container
            pair(itype,jtype)%coulomb = kCoul /= zero
            if (pair(itype,jtype)%coulomb) pair(itype,jtype)%kCoul = kCoul
            call pair(itype,jtype) % model % shifting_setup( me%Rc )
            pair(jtype,itype) = pair(itype,jtype)
          end if
        end associate
      class default
        call error( "pair model setup", "a valid pair model must be provided" )
    end select

  end subroutine set_pair_type

!===================================================================================================

  subroutine allocate_rigid_bodies( me, bodies )
    type(tData), intent(inout) :: me
    type(c_ptr), intent(in)    :: bodies

    integer :: i, j
    integer, allocatable :: seq(:), indices(:)
    integer,     pointer :: pbody(:)

    seq = [(i,i=1,me%natoms)]
    if (c_associated(bodies)) then
      call c_f_pointer( bodies, pbody, [me%natoms] )
      me%atomBody = clean_body_indices( pbody )
      me%nbodies = maxval( me%atomBody )
      me%free = pack( seq, me%atomBody == 0 )
      me%nfree = size(me%free)
      allocate( me%body(me%nbodies) )
      do i = 1, me%nbodies
        indices = pack( seq, me%atomBody == i )
        call me % body(i) % setup( indices, me%mass(indices) )
      end do
      ! Replace zeros by distinct indices:
      i = me%nbodies
      do j = 1, me%natoms
        if (me%atomBody(j) == 0) then
          i = i + 1
          me%atomBody(j) = i
        end if
      end do
    else
      me%nbodies = 0
      me%free = seq
      me%atomBody = seq
      me%nfree = me%natoms
      allocate( me%body(0) )
    end if
    me%threadFreeAtoms = (me%nfree + me%nthreads - 1)/me%nthreads
    me%threadBodies = (me%nbodies + me%nthreads - 1)/me%nthreads

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      function clean_body_indices( body ) result( index )
        integer, intent(in) :: body(:)
        integer             :: index(size(body))

        ! If body(i) <= 0 or body(i) is unique, then index(i) = 0.
        ! Otherwise, 1 <= index(i) <= number of distinc positive entries which are non-unique

        integer :: i, j, n, ibody
        logical :: not_found
        integer :: saved(size(body)), amount(size(body)), first(size(body))

        n = 0
        do i = 1, size(body)
          ibody = body(i)
          if (ibody > 0) then
            j = 0
            not_found = .true.
            do while (not_found.and.(j < n))
              j = j + 1
              not_found = saved(j) /= ibody
            end do
            if (not_found) then
              n = n + 1
              amount(n) = 1
              saved(n) = ibody
              first(n) = i
              index(i) = 0
            else
              amount(j) = amount(j) + 1
              index(i) = j
              index(first(j)) = j
            end if
          else
            index(i) = 0
          end if
        end do
        i = 0
        do j = 1, n
          if (amount(j) > 1) then
            i = i + 1
            saved(j) = i
          end if
        end do
        where (index > 0) index = saved(index)

      end function clean_body_indices
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine allocate_rigid_bodies

!===================================================================================================

  subroutine perform_initialization( me, DoF, RotDoF )
    type(tData), intent(inout) :: me
    integer,     intent(out)   :: DoF, RotDoF

    character(*), parameter :: task = "system initialization"

    integer :: i, j, layer, bodyDoF
    logical :: kspace_required

    integer, allocatable :: bodyAtom(:)

    !$omp parallel num_threads(me%nthreads)
    call update_rigid_bodies( me, omp_get_thread_num() + 1 )
    !$omp end parallel

    bodyDoF = sum(me%body%dof)
    RotDoF = bodyDoF - 3*me%nbodies
    DoF = 3*me%nfree + bodyDoF - 3

    call check_actual_interactions( me )

    ! Find out if any coulomb model requires a kspace solver:
    layer = 0
    kspace_required = .false.
    do while (.not.kspace_required.and.(layer < me%nlayers))
      layer = layer + 1
      kspace_required = me%coul(layer)%model%requires_kspace
    end do

    ! Now check if demand and supply are inconsistent:
    if (kspace_required .neqv. me%kspace_active) then
      if (kspace_required) then
        call error( task, "a kspace solver is required, but has not been defined" )
      else
        ! Deactivate kspace solver once it is not required:
        me%kspace_active = .false.
      end if
    end if

    if (me%kspace_active) then
      bodyAtom = [(me%body(i)%index,i=1,me%nbodies)]
      call me % kspace % initialize( me%nthreads, me%Rc, me%Lbox3, me%atomType, me%charge, &
                                     me%R, me%body%NP, bodyAtom, me%pair%kCoul )
      do layer = 1, me%nlayers
        call me % coul(layer) % model % kspace_setup( me%kspace%alpha )
      end do
    end if

    if (me%respa) then
      me%shortPair = me%pair
      allocate( me%shortF(3,me%natoms) )
      do i = 1, me%ntypes
        do j = 1, me%ntypes
          do layer = 1, me%nlayers
            associate ( pair => me%shortPair(i,j,layer)%model )
              pair%shifted_force = .true.
              call pair % shifting_setup( me%respaRc )
            end associate
          end do
        end do
      end do

    end if

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine update_rigid_bodies( me, thread )
        type(tData),  intent(inout) :: me
        integer,      intent(in)    :: thread

        integer :: i, j
        real(rb), allocatable :: R(:,:)

        do j = (thread - 1)*me%threadBodies + 1, min(thread*me%threadBodies, me%nbodies)
          associate(b => me%body(j))
            R = me%R(:,b%index)
            forall (i = 2:b%NP) R(:,i) = R(:,i) - me%Lbox*anint(me%invL*(R(:,i) - R(:,1)))
            call b % update( R )
            me%R(:,b%index) = R
          end associate
        end do

      end subroutine update_rigid_bodies
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine perform_initialization

!===================================================================================================

  subroutine rigid_body_forces( me, Virial )
    type(tData), intent(inout)         :: me
    real(rb),    intent(out), optional :: Virial

    real(rb) :: Wrb(me%nthreads)

    !$omp parallel num_threads(me%nthreads)
    block
      integer :: thread
      thread = omp_get_thread_num() + 1
      call compute_body_forces( thread, Wrb(thread) )
    end block
    !$omp end parallel
    if (present(Virial)) Virial = -sum(Wrb)

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
    integer,     intent(in)    :: M
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

        integer :: i, j, k, icell, ix, iy, iz, first, last, atoms_per_thread
        integer :: icoord(3)
        integer, allocatable :: head(:)

        if (make_cells) then
          first = (thread - 1)*cells_per_thread + 1
          last = min( thread*cells_per_thread, me%ncells )
          do icell = first, last
            k = icell - 1
            iz = k/MM
            j = k - iz*MM
            iy = j/M
            ix = j - iy*M
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
    real(rb) :: Rc2, Ed, Fd, r2, invR2, invR, Eij, Wij, Qi, Qj
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

  subroutine build_neighbor_lists( me, thread, Rs )
    type(tData), intent(inout) :: me
    integer,     intent(in)    :: thread
    real(rb),    intent(in)    :: Rs(3,me%natoms)

    integer  :: i, j, k, m, n, icell, jcell, npairs, itype, jtype, ibody, ipairs, middle
    integer  :: nlocal, ntotal, first, last
    real(rb) :: xRc2, xInRc2, r2
    logical  :: include(0:me%maxpairs)
    integer  :: atom(me%maxpairs), index(me%natoms)
    real(rb) :: Ri(3), Rij(3)

    integer,  allocatable :: xlist(:)
    real(rb), allocatable :: Ratom(:,:), Rvec(:,:)

    xRc2 = me%xRcSq*me%invL2
    xInRc2 = me%xRespaRcSq*me%invL2

    include = .true.
    index = 0
    npairs = 0
    associate ( neighbor => me%neighbor(thread) )
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
        Ratom = Rs(:,atom(1:ntotal))
        do k = 1, nlocal
          i = atom(k)
          first = npairs + 1
          ipairs = 0
          middle = 0
          itype = me%atomType(i)
          ibody = me%atomBody(i)
          Ri = Rs(:,i)
          xlist = index(me%excluded%item(me%excluded%first(i):me%excluded%last(i)))
          include(xlist) = .false.
          associate ( jlist => atom(k+1:ntotal),      &
                      item => neighbor%item(first:),  &
                      value => neighbor%value(first:) )
            Rvec = Ratom(:,k+1:ntotal)
            forall (m=1:size(jlist)) Rvec(:,m) = pbc(Ri - Rvec(:,m))
            do m = 1, size(jlist)
              j = jlist(m)
              jtype = me%atomType(j)
              if (include(k+m).and.(me%atomBody(j) /= ibody).and.me%interact(itype,jtype)) then
                Rij = Rvec(:,m)
                r2 = sum(Rij*Rij)
                if (r2 < xRc2) then
                  call insert_neighbor( item, value, ipairs, j, r2 )
                  if (r2 < xInRc2) middle = middle + 1
                end if
              end if
            end do
          end associate
          include(xlist) = .true.

          neighbor%first(i) = first
          neighbor%middle(i) = npairs + middle
          npairs = npairs + ipairs
          neighbor%last(i)  = npairs

        end do
        index(atom(1:ntotal)) = 0
      end do
      neighbor%count = npairs
    end associate

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      elemental function pbc( x )
        real(rb), intent(in) :: x
        real(rb)              :: pbc
        pbc = x - anint(x)
      end function pbc
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      pure subroutine insert_neighbor( item, value, npairs, j, r2 )
        integer,     intent(inout) :: item(:)
        real(rb),    intent(inout) :: value(:)
        integer,     intent(inout) :: npairs
        integer,     intent(in)    :: j
        real(rb),    intent(in)    :: r2
        integer :: k
        do k = npairs, 1, -1
          if (value(k) < r2) exit
          value(k+1) = value(k)
          item(k+1) = item(k)
        end do
        value(k+1) = r2
        item(k+1) = j
        npairs = npairs + 1
      end subroutine insert_neighbor
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine build_neighbor_lists

!===================================================================================================

  subroutine compute_pair_forces( me, thread, Rs, F, Wpair, Wcoul )
    type(tData), intent(inout) :: me
    integer,     intent(in)    :: thread
    real(rb),    intent(in)    :: Rs(3,me%natoms)
    real(rb),    intent(out)   :: F(3,me%natoms), Wpair, Wcoul

    real(rb) :: Rc2

    Rc2 = me%RcSq*me%invL2
#   include "compute.f90"

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      elemental function pbc( x )
        real(rb), intent(in) :: x
        real(rb)              :: pbc
        pbc = x - anint(x)
      end function pbc
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine compute_pair_forces

!===================================================================================================

  subroutine compute_pairs( me, thread, compute, Rs, F, Epair, Ecoul, Wpair, Wcoul )
    type(tData), intent(inout) :: me
    integer,     intent(in)    :: thread
    logical(lb), intent(in)    :: compute
    real(rb),    intent(in)    :: Rs(3,me%natoms)
    real(rb),    intent(out)   :: F(3,me%natoms), Epair, Ecoul, Wpair, Wcoul

    real(rb) :: Rc2

    Rc2 = me%RcSq*me%invL2
    if (compute) then
#     define compute
#     include "compute.f90"
#     undef compute
    else
#     include "compute.f90"
    end if

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      elemental function pbc( x )
        real(rb), intent(in) :: x
        real(rb)              :: pbc
        pbc = x - anint(x)
      end function pbc
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine compute_pairs


!===================================================================================================

  subroutine compute_short_range_forces( me, thread, Rs, F )
    type(tData), intent(inout) :: me
    integer,     intent(in)    :: thread
    real(rb),    intent(in)    :: Rs(3,me%natoms)
    real(rb),    intent(out)   :: F(3,me%natoms)

    integer  :: i, j, k, m, itype, jtype, firstAtom, lastAtom
    real(rb) :: Rc2, r2, invR2, invR, Wij, Qi, QiQj, WCij, rFc
    real(rb) :: Rij(3), Ri(3), Fi(3), Fij(3)
    logical  :: icharged, ijcharged
    real(rb), allocatable :: Rvec(:,:)

    Rc2 = me%respaRcSq*me%invL2
    F = zero
    associate ( neighbor => me%neighbor(thread) )
      firstAtom = me%cellAtom%first(me%threadCell%first(thread))
      lastAtom = me%cellAtom%last(me%threadCell%last(thread))
      do k = firstAtom, lastAtom
        i = me%cellAtom%item(k)
        itype = me%atomType(i)
        Qi = me%charge(i)
        icharged = me%charged(i)
        Ri = Rs(:,i)
        Fi = zero
        associate ( partner => me%shortPair(:,itype,me%layer), &
                    jlist => neighbor%item(neighbor%first(i):neighbor%middle(i)) )
          Rvec = Rs(:,jlist)
          forall (m=1:size(jlist)) Rvec(:,m) = pbc(Ri - Rvec(:,m))
          do m = 1, size(jlist)
            j = jlist(m)
            Rij = Rvec(:,m)
            r2 = sum(Rij*Rij)
            if (r2 < Rc2) then
              invR2 = me%invL2/r2
              invR = sqrt(invR2)
              jtype = me%atomType(j)
              ijcharged = icharged.and.me%charged(j)
              associate( pair => partner(jtype) )
                select type ( model => pair%model )
                  include "virial_compute_pair.f90"
                end select
                if (pair%model%shifted_force) then
                  rFc = pair%model%fshift/invR
                  Wij = Wij - rFc
                end if
                if (ijcharged.and.pair%coulomb) then
                  QiQj = pair%kCoul*Qi*me%charge(j)
                  WCij = QiQj*(invR - me%fshift/invR)
                  Wij = Wij + WCij
                end if
              end associate
              Fij = Wij*invR2*Rij
              Fi = Fi + Fij
              F(:,j) = F(:,j) - Fij
            end if
          end do
        end associate
        F(:,i) = F(:,i) + Fi
      end do
    end associate
    F = me%Lbox*F

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      elemental function pbc( x )
        real(rb), intent(in) :: x
        real(rb)              :: pbc
        pbc = x - anint(x)
      end function pbc
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine compute_short_range_forces

!===================================================================================================

  subroutine compute_kspace( me, compute, Rs, Elong, F )
    type(tData), intent(inout) :: me
    logical(lb), intent(in)    :: compute
    real(rb),    intent(in)    :: Rs(3,me%natoms)
    real(rb),    intent(inout) :: Elong, F(3,me%natoms)

    integer :: i, layer
    logical :: kspace_required

    kspace_required = me%coul(me%layer)%model%requires_kspace
    if (compute.or.kspace_required) call me % kspace % prepare( Rs )

    if (kspace_required) then
      call me % kspace % compute( me%layer, Elong, F )
      call me % kspace % discount_rigid_pairs( me%layer, me%Lbox3, me%R, Elong )
    end if

    if (compute) then
      associate ( E => me%threadEnergy(:,1) )
        if (kspace_required) E(me%layer) = E(me%layer) + Elong
        do i = 1, me%nlayers-1
          layer = me%other_layer(i)
          if (me%coul(layer)%model%requires_kspace) then
            call me % kspace % compute( layer, E(layer) )
            call me % kspace % discount_rigid_pairs( layer, me%Lbox3, me%R, E(layer) )
          end if
        end do
      end associate
    end if

  end subroutine compute_kspace

!===================================================================================================

  subroutine update_pairs( me, thread, Rs, DF, DE, DW, new_layer )
    type(tData), intent(in)  :: me
    integer,     intent(in)  :: thread, new_layer
    real(rb),    intent(in)  :: Rs(3,me%natoms)
    real(rb),    intent(out) :: DF(3,me%natoms), DE, DW

    integer  :: i, j, k, m, itype, jtype, firstAtom, lastAtom
    real(rb) :: Rc2, r2, invR2, invR, new_Eij, new_Wij, Eij, Wij, Qi, Qj
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

  subroutine move( me, thread, R_factor, P_factor, dt, translate, rotate, mode )
    type(tData), intent(inout) :: me
    integer,     intent(in)    :: thread
    real(rb),    intent(in)    :: R_factor, P_factor, dt
    logical(lb), intent(in)    :: translate, rotate
    integer(ib), intent(in)    :: mode

    integer :: i, j

    do i = (thread - 1)*me%threadBodies + 1, min(thread*me%threadBodies, me%nbodies)
      associate(b => me%body(i))
        if (translate) b%rcm = R_factor*b%rcm + P_factor*b%invMass*b%pcm
        if (rotate) then
          if (mode == 0) then
           call b % rotate_exact( dt )
          else
            call b % rotate_no_squish( dt, n = mode )
          end if
          forall (j=1:3) me%R(j,b%index) = b%rcm(j) + b%delta(j,:)
        end if
      end associate
    end do
    if (translate) then
      do i = (thread - 1)*me%threadFreeAtoms + 1, min(thread*me%threadFreeAtoms, me%nfree)
        j = me%free(i)
        me%R(:,j) = R_factor*me%R(:,j) + P_factor*me%P(:,j)*me%invMass(j)
      end do
    end if

  end subroutine move

!===================================================================================================

  subroutine boost( me, thread, P_factor, F_factor, translate, rotate )
    type(tData), intent(inout) :: me
    integer,     intent(in)    :: thread
    real(rb),    intent(in)    :: P_factor, F_factor
    logical(lb), intent(in)    :: translate, rotate

    integer  :: i, j
    real(rb) :: Ctau

    Ctau = two*F_factor
    if (translate) then
      do i = (thread - 1)*me%threadBodies + 1, min(thread*me%threadBodies, me%nbodies)
        associate(b => me%body(i))
          b%pcm = P_factor*b%pcm + F_factor*b%F
        end associate
      end do
      do i = (thread - 1)*me%threadFreeAtoms + 1, min(thread*me%threadFreeAtoms, me%nfree)
        j = me%free(i)
         me%P(:,j) = P_factor*me%P(:,j) + F_factor*me%F(:,j)
      end do
    end if
    if (rotate) then
      do i = (thread - 1)*me%threadBodies + 1, min(thread*me%threadBodies, me%nbodies)
        associate(b => me%body(i))
          call b%assign_momenta( P_factor*b%pi + matmul( matrix_C(b%q), Ctau*b%tau ) )
        end associate
      end do
    end if

  end subroutine boost

!===================================================================================================

end module EmDeeData
