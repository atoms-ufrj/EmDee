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

! TODO: 1) Optimize parallel performance of download in a unique omp parallel
! TODO: 2) Create indexing for having sequential body particles and free particles in arrays

module EmDeeCode

use, intrinsic :: ieee_arithmetic

use omp_lib
use lists
use math
use models
use structs
use ArBee

implicit none

character(11), parameter :: VERSION = "04 Nov 2016"

integer, parameter, private :: extra = 2000

integer, parameter, private :: ndiv = 2
integer, parameter, private :: nbcells = 62
integer, parameter, private :: nb(3,nbcells) = reshape( [ &
   1, 0, 0,    2, 0, 0,   -2, 1, 0,   -1, 1, 0,    0, 1, 0,    1, 1, 0,    2, 1, 0,   -2, 2, 0,  &
  -1, 2, 0,    0, 2, 0,    1, 2, 0,    2, 2, 0,   -2,-2, 1,   -1,-2, 1,    0,-2, 1,    1,-2, 1,  &
   2,-2, 1,   -2,-1, 1,   -1,-1, 1,    0,-1, 1,    1,-1, 1,    2,-1, 1,   -2, 0, 1,   -1, 0, 1,  &
   0, 0, 1,    1, 0, 1,    2, 0, 1,   -2, 1, 1,   -1, 1, 1,    0, 1, 1,    1, 1, 1,    2, 1, 1,  &
  -2, 2, 1,   -1, 2, 1,    0, 2, 1,    1, 2, 1,    2, 2, 1,   -2,-2, 2,   -1,-2, 2,    0,-2, 2,  &
   1,-2, 2,    2,-2, 2,   -2,-1, 2,   -1,-1, 2,    0,-1, 2,    1,-1, 2,    2,-1, 2,   -2, 0, 2,  &
  -1, 0, 2,    0, 0, 2,    1, 0, 2,    2, 0, 2,   -2, 1, 2,   -1, 1, 2,    0, 1, 2,    1, 1, 2,  &
   2, 1, 2,   -2, 2, 2,   -1, 2, 2,    0, 2, 2,    1, 2, 2,    2, 2, 2 ], shape(nb) )

type, bind(C) :: tEmDee
  integer(ib) :: builds         ! Number of neighbor list builds
  real(rb)    :: pairTime       ! Time taken in force calculations
  real(rb)    :: totalTime      ! Total time since initialization
  real(rb)    :: Potential      ! Total potential energy of the system
  real(rb)    :: Kinetic        ! Total kinetic energy of the system
  real(rb)    :: Rotational     ! Rotational kinetic energy of the system
  real(rb)    :: Virial         ! Total internal virial of the system
  integer(ib) :: DOF            ! Total number of degrees of freedom
  integer(ib) :: RDOF           ! Number of rotational degrees of freedom
  integer(ib) :: rotationMode   ! Algorithm used for free rotation of rigid bodies
  type(c_ptr) :: Data           ! Pointer to system data
end type tEmDee

type, private :: tCell
  integer :: neighbor(nbcells)
end type tCell

type, private :: tData

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

  logical :: coulomb = .false.            ! Flag for checking if coulombic interactions occur
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
  real(rb), allocatable :: charge(:,:)    ! Electric charges of all atoms
  real(rb), allocatable :: mass(:)        ! Masses of all atoms
  real(rb), allocatable :: invMass(:)     ! Inverses of atoms masses
  real(rb), allocatable :: R0(:,:)        ! Position of each atom at latest neighbor list building

  type(c_ptr), allocatable :: model(:,:,:) ! Model of each type of atom pair
  type(tCell), allocatable :: cell(:)      ! Array containing all neighbor cells of each cell
  type(tBody), allocatable :: body(:)      ! Pointer to the rigid bodies present in the system
  type(tList), allocatable :: neighbor(:)  ! Pointer to neighbor lists

end type tData

private :: rigid_body_forces, maximum_approach_sq, distribute_atoms, find_pairs_and_compute, &
           compute_pairs, compute_bonds, compute_angles, compute_dihedrals, compute_group_energy

contains

!===================================================================================================
!                                L I B R A R Y   P R O C E D U R E S
!===================================================================================================

  function EmDee_system( threads, layers, rc, skin, N, types, masses ) bind(C,name="EmDee_system")
    integer(ib), value :: threads, layers, N
    real(rb),    value :: rc, skin
    type(c_ptr), value :: types, masses
    type(tEmDee)       :: EmDee_system

    integer :: i
    integer,     pointer :: ptype(:)
    real(rb),    pointer :: pmass(:)
    type(tData), pointer :: me

    write(*,'("EmDee (version: ",A11,")")') VERSION

    ! Allocate data structure:
    allocate( me )

    ! Set up fixed entities:
    me%nthreads = threads
    me%nlayers = layers
    me%Rc = rc
    me%RcSq = rc*rc
    me%xRc = rc + skin
    me%xRcSq = me%xRc**2
    me%skinSq = skin*skin
    me%natoms = N
    me%fshift = one/me%RcSq
    me%eshift = -two/me%Rc

    ! Set up atom types:
    if (c_associated(types)) then
      call c_f_pointer( types, ptype, [N] )
      if (minval(ptype) /= 1) stop "ERROR: wrong specification of atom types."
      me%ntypes = maxval(ptype)
      allocate( me%atomType(N), source = ptype )
    else
      me%ntypes = 1
      allocate( me%atomType(N), source = 1 )
    end if

    ! Set up atomic masses:
    if (c_associated(masses)) then
      call c_f_pointer( masses, pmass, [me%ntypes] )
      allocate( me%mass(N), source = pmass(ptype) )
      allocate( me%invMass(N), source = one/pmass(ptype) )
      me%totalMass = sum(pmass(ptype))
    else
      allocate( me%mass(N), source = one )
      allocate( me%invMass(N), source = one )
      me%totalMass = real(N,rb)
    end if

    ! Initialize counters and other mutable entities:
    me%startTime = omp_get_wtime()
    allocate( me%R(3,N), me%P(3,N), me%F(3,N), me%R0(3,N), source = zero )
    allocate( me%charge(N,layers), source = zero )
    allocate( me%cell(0) )
    allocate( me%atomCell(N) )

    ! Allocate variables associated to rigid bodies:
    allocate( me%body(me%nbodies) )
    me%nfree = N
    allocate( me%free(N), source = [(i,i=1,N)] )
    me%threadAtoms = (N + threads - 1)/threads

    ! Allocate memory for list of atoms per cell:
    call me % cellAtom % allocate( N, 0 )

    ! Allocate memory for neighbor lists:
    allocate( me%neighbor(threads) )
    call me % neighbor % allocate( extra, N )

    ! Allocate memory for the list of pairs excluded from the neighbor lists:
    call me % excluded % allocate( extra, N )

    ! Allocate memory for pair models:
    allocate( me%model(me%ntypes,me%ntypes,me%nlayers) )
    me%model = c_null_ptr

    ! Set up mutable entities:
    EmDee_system % builds = 0
    EmDee_system % pairTime = zero
    EmDee_system % totalTime = zero
    EmDee_system % Potential = zero
    EmDee_system % Kinetic = zero
    EmDee_system % Rotational = zero
    EmDee_system % DOF = 3*(N - 1)
    EmDee_system % RDOF = 0
    EmDee_system % rotationMode = 0
    EmDee_system % data = c_loc(me)

  end function EmDee_system

!===================================================================================================

  subroutine EmDee_switch_model_layer( md, layer ) bind(C,name="EmDee_set_layer")
    type(tEmDee), intent(inout) :: md
    integer(ib),  value         :: layer
    type(tData), pointer :: me
    call c_f_pointer( md%data, me )
    if ((layer < 1).or.(layer > me%nlayers)) stop "ERROR in model layer change: out of range"
    me%layer = layer
    if (me%initialized) call compute_forces( md )
  end subroutine EmDee_switch_model_layer

!===================================================================================================

  subroutine EmDee_set_charges( md, charges ) bind(C,name="EmDee_set_charges")
    type(tEmDee), value :: md
    type(c_ptr),  value :: charges

    integer :: i, j
    type(tData),  pointer :: me
    type(tModel), pointer :: model
    real(rb),     pointer :: Q(:)

    call c_f_pointer( md%data, me )
    if (.not.me%coulomb) then
      me%coulomb = .true.
      do i = 1, me%ntypes
        do j = 1, me%ntypes
          call c_f_pointer( me%model(i,j,me%layer), model )
          if (associated(model)) model%id = mCOULOMB + mod(model%id,mCOULOMB)
        end do
      end do
    end if
    call c_f_pointer( charges, Q, [me%natoms] )
    me%charge(:,me%layer) = Q
    if (me%initialized) call compute_forces( md )

  end subroutine EmDee_set_charges

!===================================================================================================

  subroutine EmDee_set_pair_type( md, itype, jtype, model ) bind(C,name="EmDee_set_pair_type")
    type(tEmDee), value :: md
    integer(ib),  value :: itype, jtype
    type(c_ptr),  value :: model

    integer :: k
    logical :: keep
    type(tData),  pointer :: me
    type(tModel), pointer :: ikModel

    call c_f_pointer( md%data, me )
    if (me%initialized) stop "ERROR: cannot set pair type after coordinates have been provided"
    if (itype == jtype) then
      call associate_model( itype, itype, model )
      do k = 1, me%ntypes
        if (k /= itype) then
          call c_f_pointer( me%model(itype,k,me%layer), ikModel )
          keep = associated(ikModel)
          if (keep) keep = ikModel%external
          if (.not.keep) call replace( itype, k, cross_pair( model, me%model(k,k,me%layer) ) )
        end if
      end do
    else
      call replace( itype, jtype, model )
    end if

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine associate_model( i, j, model )
        integer, intent(in) :: i, j
        type(c_ptr), intent(in) :: model
        type(tModel), pointer :: pmodel
        if (c_associated(model)) then
          if (me%coulomb) then
            call c_f_pointer( model, pmodel )
            pmodel%id = mCOULOMB + mod(pmodel%id,mCOULOMB)
          end if
          me%model(i,j,me%layer) = model
        else
          me%model(i,j,me%layer) = c_null_ptr
        end if
      end subroutine associate_model
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine replace( i, j, model )
        integer, intent(in) :: i, j
        type(c_ptr), intent(in) :: model
        type(tModel), pointer :: ij
        if (c_associated(me%model(i,j,me%layer))) then
          call c_f_pointer( me%model(i,j,me%layer), ij )
          if (.not.ij%external) deallocate( ij )
        end if
        call associate_model( i, j, model )
        call associate_model( j, i, model )
      end subroutine replace
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine EmDee_set_pair_type

!===================================================================================================

  subroutine EmDee_ignore_pair( md, i, j ) bind(C,name="EmDee_ignore_pair")
    type(tEmDee), value :: md
    integer(ib),  value :: i, j

    integer :: n
    type(tData), pointer :: me

    call c_f_pointer( md%data, me )
    if ((i > 0).and.(i <= me%natoms).and.(j > 0).and.(j <= me%natoms).and.(i /= j)) then
      associate (excluded => me%excluded)
        n = excluded%count
        if (n == excluded%nitems) call excluded % resize( n + extra )
        call add_item( excluded, i, j )
        call add_item( excluded, j, i )
        excluded%count = n
      end associate
    end if

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine add_item( excluded, i, j )
        type(tList), intent(inout) :: excluded
        integer,     intent(in)    :: i, j
        integer :: start, end
        start = excluded%first(i)
        end = excluded%last(i)
#define core1 excluded%item(end+2:n+1) = excluded%item(end+1:n)
#define core2 excluded%item(end+1) = j
#define core core1; core2
!ORIGINAL VERSION SUBJECT TO A FALSE ALARM ERROR BY BOUNDS-CHECK
!       if ((end < start).or.(j > excluded%item(end))) then
!         excluded%item(end+2:n+1) = excluded%item(end+1:n)
!         excluded%item(end+1) = j
! ALTERNATIVE USING -CPP #DEFINE AND #UNDEF
        if (end < start) then
          core
        elseif (j > excluded%item(end)) then
          core
!END ALTERNATIVE
        else
          do while (j > excluded%item(start))
            start = start + 1
          end do
          if (j == excluded%item(start)) return
          excluded%item(start+1:n+1) = excluded%item(start:n)
          excluded%item(start) = j
          start = start + 1
        end if
#undef core
#undef core2
#undef core1
        excluded%first(i+1:) = excluded%first(i+1:) + 1
        excluded%last(i:) = excluded%last(i:) + 1
        n = n + 1
      end subroutine add_item
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine EmDee_ignore_pair

!===================================================================================================

  subroutine EmDee_add_bond( md, i, j, model ) bind(C,name="EmDee_add_bond")
    type(tEmDee), value :: md
    integer(ib),  value :: i, j
    type(c_ptr),  value :: model

    type(tData), pointer :: me

    call c_f_pointer( md%data, me )
    call me % bonds % add( i, j, 0, 0, model )
    call EmDee_ignore_pair( md, i, j )

  end subroutine EmDee_add_bond

!===================================================================================================

  subroutine EmDee_add_angle( md, i, j, k, model ) bind(C,name="EmDee_add_angle")
    type(tEmDee), value :: md
    integer(ib),  value :: i, j, k
    type(c_ptr),  value :: model

    type(tData), pointer :: me

    call c_f_pointer( md%data, me )
    call me % angles % add( i, j, k, 0, model )
    call EmDee_ignore_pair( md, i, j )
    call EmDee_ignore_pair( md, i, k )
    call EmDee_ignore_pair( md, j, k )

  end subroutine EmDee_add_angle

!===================================================================================================

  subroutine EmDee_add_dihedral( md, i, j, k, l, model ) bind(C,name="EmDee_add_dihedral")
    type(tEmDee), value :: md
    integer(ib),  value :: i, j, k, l
    type(c_ptr),  value :: model

    type(tData), pointer :: me

    call c_f_pointer( md%data, me )
    call me % dihedrals % add( i, j, k, l, model )
    call EmDee_ignore_pair( md, i, j )
    call EmDee_ignore_pair( md, i, k )
    call EmDee_ignore_pair( md, i, l )
    call EmDee_ignore_pair( md, j, k )
    call EmDee_ignore_pair( md, j, l )
    call EmDee_ignore_pair( md, k, l )

  end subroutine EmDee_add_dihedral

!===================================================================================================

  subroutine EmDee_add_rigid_body( md, N, indexes ) bind(C,name="EmDee_add_rigid_body")
    type(tEmDee), value :: md
    type(c_ptr),  value :: indexes
    integer(ib),  value :: N

    integer, parameter :: extraBodies = 100

    integer :: i, j
    logical,  allocatable :: isFree(:)
    real(rb), allocatable :: Rn(:,:)
    integer,      pointer :: atom(:)
    type(tData),  pointer :: me

    type(tBody), allocatable :: body(:)

    call c_f_pointer( indexes, atom, [N] )
    call c_f_pointer( md%data, me )

    allocate( isFree(me%natoms) )
    isFree = .false.
    isFree(me%free(1:me%nfree)) = .true.
    isFree(atom) = .false.
    me%nfree = me%nfree - N
    md%DOF = md%DOF - 3*N
    if (count(isFree) /= me%nfree) stop "Error adding rigid body: only free atoms are allowed."
    me%free(1:me%nfree) = pack([(i,i=1,me%natoms)],isFree)
    me%threadAtoms = (me%nfree + me%nthreads - 1)/me%nthreads

    if (me%nbodies == size(me%body)) then
      allocate( body(me%nbodies + extraBodies) )
      body(1:me%nbodies) = me%body
      deallocate( me%body )
      call move_alloc( body, me%body )
    end if
    me%nbodies = me%nbodies + 1
    me%threadBodies = (me%nbodies + me%nthreads - 1)/me%nthreads
    associate(b => me%body(me%nbodies))
      call b % setup( atom, me%mass(atom) )
      if (me%initialized) then
        Rn = me%R(:,atom)
        forall (j=2:b%NP) Rn(:,j) = Rn(:,j) - me%Lbox*anint((Rn(:,j) - Rn(:,1))*me%invL)
        call b % update( Rn )
        me%R(:,b%index) = Rn
      end if
      md%DOF = md%DOF + b%dof
      md%RDOF = md%RDOF + b%dof - 3
    end associate
    do i = 1, N-1
      do j = i+1, N
        call EmDee_ignore_pair( md, atom(i), atom(j) )
      end do
    end do

  end subroutine EmDee_add_rigid_body

!===================================================================================================

  subroutine EmDee_upload( md, Lbox, coords, momenta, forces ) bind(C,name="EmDee_upload")
    type(tEmDee), intent(inout) :: md
    type(c_ptr),  value         :: Lbox, coords, momenta, forces

    real(rb) :: twoKEt, twoKEr
    real(rb), pointer :: L
    type(tData), pointer :: me

    call c_f_pointer( md%data, me )

    if (.not.me%initialized) then
      me%initialized = c_associated(Lbox) .and. c_associated(coords)
      if (.not.me%initialized) then
        stop "ERROR in EmDee_upload: box side length and atomic coordinates are required."
      end if
    end if

    if (c_associated(Lbox)) then
      call c_f_pointer( Lbox, L )
      me%Lbox = L
      me%invL = one/L
      me%invL2 = me%invL**2
    end if

    !$omp parallel num_threads(me%nthreads) reduction(+:twoKEt,twoKEr)
    block
      integer :: thread
      thread = omp_get_thread_num() + 1
      if (c_associated(coords)) call assign_coordinates( thread )
      if (c_associated(momenta)) call assign_momenta( thread, twoKEt, twoKEr )
      if (c_associated(forces)) call assign_forces( thread )
    end block
    !$omp end parallel

    if (c_associated(coords).and.(.not.c_associated(forces))) call compute_forces( md )
    if (c_associated(momenta)) then
      md%Rotational = half*twoKEr
      md%Kinetic = half*twoKEt + md%Rotational
    end if

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine assign_coordinates( thread )
        integer, intent(in) :: thread
        integer :: i, j
        real(rb), allocatable :: R(:,:)
        real(rb), pointer :: Rext(:,:)
        call c_f_pointer( coords, Rext, [3,me%natoms] )
        do j = (thread - 1)*me%threadAtoms + 1, min(thread*me%threadAtoms, me%nfree)
          i = me%free(j)
          me%R(:,i) = Rext(:,i)
        end do
        do i = (thread - 1)*me%threadBodies + 1, min(thread*me%threadBodies, me%nbodies)
          associate(b => me%body(i))
            R = Rext(:,b%index)
            forall (j=2:b%NP) R(:,j) = R(:,j) - me%Lbox*anint((R(:,j) - R(:,1))*me%invL)
            call b % update( R )
            me%R(:,b%index) = R
          end associate
        end do
      end subroutine assign_coordinates
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine assign_momenta( thread, twoKEt, twoKEr )
        integer,  intent(in)  :: thread
        real(rb), intent(out) :: twoKEt, twoKEr
        integer :: i, j
        real(rb) :: L(3), Pj(3)
        real(rb), pointer :: Pext(:,:)
        twoKEt = zero
        twoKEr = zero
        call c_f_pointer( momenta, Pext, [3,me%natoms] )
        do j = (thread - 1)*me%threadAtoms + 1, min(thread*me%threadAtoms, me%nfree)
          i = me%free(j)
          me%P(:,i) = Pext(:,i)
          twoKEt = twoKEt + me%invMass(i)*sum(Pext(:,i)**2)
        end do
        do i = (thread - 1)*me%threadBodies + 1, min(thread*me%threadBodies, me%nbodies)
          associate(b => me%body(i))
            b%pcm = zero
            L = zero
            do j = 1, b%NP
              Pj = Pext(:,b%index(j))
              b%pcm = b%pcm + Pj
              L = L + cross_product( b%delta(:,j), Pj )
            end do
            twoKEt = twoKEt + b%invMass*sum(b%pcm**2)
            b%pi = matmul( matrix_C(b%q), two*L )
            b%omega= half*b%invMoI*matmul( matrix_Bt(b%q), b%pi )
            twoKEr = twoKEr + sum(b%MoI*b%omega**2)
          end associate
        end do
      end subroutine assign_momenta
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine assign_forces( thread )
        integer, intent(in) :: thread
        integer :: i, j
        real(rb), pointer :: Fext(:,:)
        call c_f_pointer( forces, Fext, [3,me%natoms] )
        do j = (thread - 1)*me%threadAtoms + 1, min(thread*me%threadAtoms, me%nfree)
          i = me%free(j)
          me%F(:,i) = Fext(:,i)
        end do
        do i = (thread - 1)*me%threadBodies + 1, min(thread*me%threadBodies, me%nbodies)
          call me % body(i) % force_torque_virial( Fext )
        end do
      end subroutine assign_forces
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine EmDee_upload

!===================================================================================================

  subroutine EmDee_download( md, Lbox, coords, momenta, forces ) bind(C,name="EmDee_download")
    type(tEmDee), value :: md
    type(c_ptr),  value :: Lbox, coords, momenta, forces

    real(rb),    pointer :: L, Pext(:,:), Rext(:,:), Fext(:,:)
    type(tData), pointer :: me

    call c_f_pointer( md%data, me )

    if (c_associated(Lbox)) then
      call c_f_pointer( Lbox, L )
      L = me%Lbox
    end if

    if (c_associated(coords)) then
      call c_f_pointer( coords, Rext, [3,me%natoms] )
      Rext = me%R
    end if

    if (c_associated(forces)) then
      call c_f_pointer( forces, Fext, [3,me%natoms] )
      Fext = me%F
    end if

    if (c_associated(momenta)) then
      call c_f_pointer( momenta, Pext, [3,me%natoms] )
      !$omp parallel num_threads(me%nthreads)
      call get_momenta( omp_get_thread_num() + 1 )
      !$omp end parallel
    end if

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine get_momenta( thread )
        integer, intent(in) :: thread
        integer :: i
        forall (i = (thread - 1)*me%threadAtoms + 1 : min(thread*me%threadAtoms, me%nfree))
          Pext(:,me%free(i)) = me%P(:,me%free(i))
        end forall
        forall(i = (thread - 1)*me%threadBodies + 1 : min(thread*me%threadBodies,me%nbodies))
          Pext(:,me%body(i)%index) = me%body(i) % particle_momenta()
        end forall
      end subroutine get_momenta
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine EmDee_download

!===================================================================================================

  subroutine EmDee_random_momenta( md, kT, adjust, seed ) bind(C,name="EmDee_random_momenta")
    type(tEmDee), intent(inout) :: md
    real(rb),     value         :: kT
    integer(ib),  value         :: adjust, seed

    integer  :: i, j
    real(rb) :: twoKEt, TwoKEr
    type(tData), pointer :: me

    call c_f_pointer( md%data, me )
    if (me%random%seeding_required) call me % random % setup( seed )
    twoKEt = zero
    TwoKEr = zero
    associate (rng => me%random)
      if (me%nbodies /= 0) then
        if (.not.me%initialized) stop "ERROR in random momenta: coordinates not defined."
        do i = 1, me%nbodies
          associate (b => me%body(i))
            b%pcm = sqrt(b%mass*kT)*[rng%normal(), rng%normal(), rng%normal()]
            call b%assign_momenta( sqrt(b%invMoI*kT)*[rng%normal(), rng%normal(), rng%normal()] )
            twoKEt = twoKEt + b%invMass*sum(b%pcm*b%pcm)
            TwoKEr = TwoKEr + sum(b%MoI*b%omega**2)
          end associate
        end do
      end if
      do j = 1, me%nfree
        i = me%free(j)
        me%P(:,i) = sqrt(me%mass(i)*kT)*[rng%normal(), rng%normal(), rng%normal()]
        twoKEt = twoKEt + sum(me%P(:,i)**2)/me%mass(i)
      end do
    end associate
    if (adjust == 1) call adjust_momenta
    md%Rotational = half*TwoKEr
    md%Kinetic = half*(twoKEt + TwoKEr)

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine adjust_momenta
        integer  :: i
        real(rb) :: vcm(3), factor
        associate (free => me%free(1:me%nfree), body => me%body(1:me%nbodies))
          forall (i=1:3) vcm(i) = (sum(me%P(i,free)) + sum(body(1:me%nbodies)%pcm(i)))/me%totalMass
          forall (i=1:me%nfree) me%P(:,free(i)) = me%P(:,free(i)) - me%mass(free(i))*vcm
          forall (i=1:me%nbodies) body(i)%pcm = body(i)%pcm - body(i)%mass*vcm
          twoKEt = sum([(sum(me%P(:,free(i))**2)*me%invMass(free(i)),i=1,me%nfree)]) + &
                   sum([(sum(body(i)%pcm**2)*body(i)%invMass,i=1,me%nbodies)])
          factor = sqrt((3*me%nfree + sum(body%dof) - 3)*kT/(twoKEt + TwoKEr))
          me%P(:,free) = factor*me%P(:,free)
          do i = 1, me%nbodies
            associate( b => body(i) )
              b%pcm = factor*b%pcm
              call b%assign_momenta( factor*b%omega )
            end associate
          end do
        end associate
        twoKEt = factor*factor*twoKEt
        TwoKEr = factor*factor*TwoKEr
      end subroutine adjust_momenta
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine EmDee_random_momenta

!===================================================================================================

  subroutine EmDee_save_state( md, rigid )
    type(tEmDee), intent(inout) :: md
    integer(ib),  intent(in)    :: rigid
    if (rigid /= 0) then
    else
    end if
  end subroutine EmDee_save_state

!===================================================================================================

  subroutine EmDee_restore_state( md )
    type(tEmDee), intent(inout) :: md
  end subroutine EmDee_restore_state

!===================================================================================================

  subroutine EmDee_boost( md, lambda, alpha, dt, translation, rotation ) bind(C,name="EmDee_boost")
    type(tEmDee), intent(inout) :: md
    real(rb),     value         :: lambda, alpha, dt
    integer(ib),  value         :: translation, rotation

    real(rb) :: CP, CF, Ctau, twoKEt, twoKEr, KEt
    logical  :: tflag, rflag
    type(tData), pointer :: me

    call c_f_pointer( md%data, me )

    CF = phi(alpha*dt)*dt
    CP = one - alpha*CF
    CF = lambda*CF
    Ctau = two*CF

    tflag = translation /= 0
    rflag = rotation /= 0
    twoKEt = zero
    twoKEr = zero
    !$omp parallel num_threads(me%nthreads) reduction(+:twoKEt,twoKEr)
    call boost( omp_get_thread_num() + 1, twoKEt, twoKEr )
    !$omp end parallel
    if (tflag) then
      KEt = half*twoKEt
    else
      KEt = md%Kinetic - md%Rotational
    end if
    if (rflag) md%Rotational = half*twoKEr
    md%Kinetic = KEt + md%Rotational

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine boost( thread, twoKEt, twoKEr )
        integer,  intent(in)    :: thread
        real(rb), intent(inout) :: twoKEt, twoKEr
        integer  :: i, j
        do i = (thread - 1)*me%threadBodies + 1, min(thread*me%threadBodies, me%nbodies)
          associate(b => me%body(i))
            if (tflag) then
              b%pcm = CP*b%pcm + CF*b%F
              twoKEt = twoKEt + b%invMass*sum(b%pcm*b%pcm)
            end if
            if (rflag) then
              call b%assign_momenta( CP*b%pi + matmul( matrix_C(b%q), Ctau*b%tau ) )
              twoKEr = twoKEr + sum(b%MoI*b%omega*b%Omega)
            end if
          end associate
        end do
        if (tflag) then
          do i = (thread - 1)*me%threadAtoms + 1, min(thread*me%threadAtoms, me%nfree)
            j = me%free(i)
            me%P(:,j) = CP*me%P(:,j) + CF*me%F(:,j)
            twoKEt = twoKEt + me%invMass(j)*sum(me%P(:,j)**2)
          end do
        end if
      end subroutine boost
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine EmDee_boost

!===================================================================================================

  subroutine EmDee_move( md, lambda, alpha, dt ) bind(C,name="EmDee_move")
    type(tEmDee), intent(inout) :: md
    real(rb),     value         :: lambda, alpha, dt

    real(rb) :: cR, cP
    type(tData), pointer :: me

    call c_f_pointer( md%data, me )

    if (alpha /= zero) then
      cP = phi(alpha*dt)*dt
      cR = one - alpha*cP
      me%Lbox = cR*me%Lbox
      me%InvL = one/me%Lbox
      me%invL2 = me%invL*me%invL
    else
      cP = dt
      cR = one
    end if
    cP = lambda*cP

    !$omp parallel num_threads(me%nthreads)
    call move( omp_get_thread_num() + 1, cP, cR )
    !$omp end parallel

    call compute_forces( md )

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine move( thread, cP, cR )
        integer,  intent(in) :: thread
        real(rb), intent(in) :: cP, cR
        integer :: i, j
        do i = (thread - 1)*me%threadBodies + 1, min(thread*me%threadBodies, me%nbodies)
          associate(b => me%body(i))
            b%rcm = cR*b%rcm + cP*b%invMass*b%pcm
            if (md%rotationMode == 0) then
              call b % rotate_exact( dt )
            else
              call b % rotate_approx( dt, n = md%rotationMode )
            end if
            forall (j=1:3) me%R(j,b%index) = b%rcm(j) + b%delta(j,:)
          end associate
        end do
        do i = (thread - 1)*me%threadAtoms + 1, min(thread*me%threadAtoms, me%nfree)
          j = me%free(i)
          me%R(:,j) = cR*me%R(:,j) + cP*me%P(:,j)*me%invMass(j)
        end do
      end subroutine move
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine EmDee_move

!===================================================================================================

  subroutine EmDee_group_energy( md, na, atoms, energies ) bind(C,name="EmDee_group_energy")
    type(tEmDee), value :: md
    integer(ib),  value :: na
    type(c_ptr),  value :: atoms, energies

    real(rb), allocatable :: energy(:)
    integer,     pointer :: atom(:)
    real(rb),    pointer :: Eext(:)
    type(tData), pointer :: me

    call c_f_pointer( md%data, me )
    call c_f_pointer( atoms, atom, [na] )
    allocate( energy(me%nlayers) )

    !$omp parallel num_threads(me%nthreads) reduction(+:energy)
    call compute_group_energy( me, omp_get_thread_num() + 1, na, atom, energy )
    !$omp end parallel

    call c_f_pointer( energies, Eext, [me%nlayers] )
    Eext = energy

  end subroutine EmDee_group_energy

!===================================================================================================
!                              A U X I L I A R Y   P R O C E D U R E S
!===================================================================================================

  subroutine compute_forces( md )
    type(tEmDee), intent(inout) :: md

    integer  :: M
    real(rb) :: time, E, W
    logical  :: buildList
    real(rb), allocatable :: Rs(:,:), Fs(:,:,:)
    type(tData), pointer :: me

    call c_f_pointer( md%data, me )
    md%pairTime = md%pairTime - omp_get_wtime()

    allocate( Rs(3,me%natoms), Fs(3,me%natoms,me%nthreads) )
    Rs = me%invL*me%R
    Fs = zero
    E = zero
    W = zero

    buildList = maximum_approach_sq( me%natoms, me%R - me%R0 ) > me%skinSq
    if (buildList) then
      M = floor(ndiv*me%Lbox/me%xRc)
      call distribute_atoms( me, max(M,2*ndiv+1), Rs )
      me%R0 = me%R
      md%builds = md%builds + 1
    endif

    !$omp parallel num_threads(me%nthreads) reduction(+:E,W)
    block
      integer :: thread
      thread = omp_get_thread_num() + 1
      associate( F => Fs(:,:,thread) )
        if (buildList) then
          call find_pairs_and_compute( me, thread, Rs, F, E, W )
        else
          call compute_pairs( me, thread, Rs, F, E, W )
        end if
        if (me%bonds%exist) call compute_bonds( me, thread, Rs, F, E, W )
        if (me%angles%exist) call compute_angles( me, thread, Rs, F, E, W )
        if (me%dihedrals%exist) call compute_dihedrals( me, thread, Rs, F, E, W )
      end associate
    end block
    !$omp end parallel

    me%F = me%Lbox*sum(Fs,3)
    md%Potential = E
    md%Virial = third*W
    if (me%nbodies /= 0) call rigid_body_forces( me, md%Virial )

    time = omp_get_wtime()
    md%pairTime = md%pairTime + time
    md%totalTime = time - me%startTime

  end subroutine compute_forces

!===================================================================================================

  subroutine compute_forces_old( md )
    type(tEmDee), intent(inout) :: md

    integer  :: M
    real(rb) :: Potential, Virial, time
    logical  :: buildList
    real(rb), allocatable :: Rs(:,:), Fs(:,:)
    type(tData),  pointer :: me

    call c_f_pointer( md%data, me )
    md%pairTime = md%pairTime - omp_get_wtime()

    allocate( Rs(3,me%natoms), Fs(3,me%natoms) )
    Rs = me%invL*me%R
    Fs = zero
    Potential = zero
    Virial = zero

    buildList = maximum_approach_sq( me%natoms, me%R - me%R0 ) > me%skinSq
    if (buildList) then
      M = floor(ndiv*me%Lbox/me%xRc)
      call distribute_atoms( me, max(M,2*ndiv+1), Rs )
      me%R0 = me%R
      md%builds = md%builds + 1
    endif

    !$omp parallel num_threads(me%nthreads) reduction(+:Fs,Potential,Virial)
    block
      integer :: thread
      thread = omp_get_thread_num() + 1
      if (buildList) then
        call find_pairs_and_compute( me, thread, Rs, Fs, Potential, Virial )
      else
        call compute_pairs( me, thread, Rs, Fs, Potential, Virial )
      end if
      if (me%bonds%exist) call compute_bonds( me, thread, Rs, Fs, Potential, Virial )
      if (me%angles%exist) call compute_angles( me, thread, Rs, Fs, Potential, Virial )
      if (me%dihedrals%exist) call compute_dihedrals( me, thread, Rs, Fs, Potential, Virial )
    end block
    !$omp end parallel

    me%F = me%Lbox*Fs
    md%Potential = Potential
    md%Virial = third*Virial
    if (me%nbodies /= 0) call rigid_body_forces( me, md%Virial )

    time = omp_get_wtime()
    md%pairTime = md%pairTime + time
    md%totalTime = time - me%startTime

  end subroutine compute_forces_old

!===================================================================================================

  subroutine rigid_body_forces( me, Virial )
    type(tData), intent(inout) :: me
    real(rb),    intent(inout) :: Virial

    real(rb) :: Wrb

    !$omp parallel num_threads(me%nthreads) reduction(+:Wrb)
    call compute_body_forces( omp_get_thread_num() + 1, Wrb )
    !$omp end parallel
    Virial = Virial - third*Wrb

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
    integer, intent(in) :: N
    real(rb),    intent(in) :: delta(3,N)

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

    integer  :: i, j, m, nbonds
    real(rb) :: d, E, mdEdr
    real(rb) :: Rij(3), Fij(3)
    type(tModel), pointer :: model

    nbonds = (me%bonds%number + me%nthreads - 1)/me%nthreads
    do m = (threadId - 1)*nbonds + 1, min( me%bonds%number, threadId*nbonds )
      associate(bond => me%bonds%item(m))
        i = bond%i
        j = bond%j
        model => bond%model
      end associate
      Rij = R(:,i) - R(:,j)
      Rij = Rij - anint(Rij)
      d = me%Lbox*sqrt(sum(Rij*Rij))
      call compute_bond
      Potential = Potential + E
      Virial = Virial + mdEdr*d
      Fij = mdEdr*Rij/d
      F(:,i) = F(:,i) + Fij
      F(:,j) = F(:,j) - Fij
    end do

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      include "compute_bond.f90"
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
    type(tModel), pointer :: model

    nangles = (me%angles%number + me%nthreads - 1)/me%nthreads
    do m = (threadId - 1)*nangles + 1, min( me%angles%number, threadId*nangles )
      associate (angle => me%angles%item(m))
        i = angle%i
        j = angle%j
        k = angle%k
        model => angle%model
      end associate
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
      call compute_angle()
      Fa = Fa/(me%Lbox*axb)
      Fi = Fa*(b - (ab/aa)*a)
      Fk = Fa*(a - (ab/bb)*b)
      F(:,i) = F(:,i) + Fi
      F(:,k) = F(:,k) + Fk
      F(:,j) = F(:,j) - (Fi + Fk)
      Potential = Potential + Ea
      Virial = Virial + me%Lbox*sum(Fi*a + Fk*b)
    end do

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      include "compute_angle.f90"
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine compute_angles

!===================================================================================================

  subroutine compute_dihedrals( me, threadId, R, F, Potential, Virial )
    type(tData), intent(inout) :: me
    integer,     intent(in)    :: threadId
    real(rb),    intent(in)    :: R(3,me%natoms)
    real(rb),    intent(inout) :: F(3,me%natoms), Potential, Virial

    integer  :: m, ndihedrals, i, j
    real(rb) :: Rc2, Ed, Fd, r2, invR2, Eij, Wij, icharge
    real(rb) :: Rj(3), Rk(3), Fi(3), Fk(3), Fl(3), Fij(3)
    real(rb) :: normRkj, normX, a, b, phi
    real(rb) :: rij(3), rkj(3), rlk(3), x(3), y(3), z(3), u(3), v(3), w(3)
    type(tModel), pointer :: model

    Rc2 = me%RcSq*me%invL2
    ndihedrals = (me%dihedrals%number + me%nthreads - 1)/me%nthreads
    do m = (threadId - 1)*ndihedrals + 1, min( me%dihedrals%number, threadId*ndihedrals )
      associate(d => me%dihedrals%item(m))
        Rj = R(:,d%j)
        Rk = R(:,d%k)
        rij = R(:,d%i) - Rj
        rkj = Rk - Rj
        rlk = R(:,d%l) - Rk
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
        model => d%model
        call compute_dihedral()
        Fd = Fd/(me%Lbox*(a*a + b*b))
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
        Potential = Potential + Ed
        Virial = Virial + me%Lbox*sum(Fi*rij + Fk*rkj + Fl*(rlk + rkj))
        if (model%p1 /= zero) then
          i = d%i
          j = d%l
          rij = rij + rlk - rkj
          r2 = sum(rij*rij)
          if (r2 < me%RcSq) then
            invR2 = me%invL2/r2
            call c_f_pointer( me%model(me%atomType(i),me%atomType(j),me%layer), model )
            icharge = me%charge(i,me%layer)
            call compute_pair()
            Eij = model%p1*Eij
            Wij = model%p1*Wij
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
      include "compute_dihedral.f90"
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      include "compute_pair.f90"
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine compute_dihedrals

!===================================================================================================

  subroutine find_pairs_and_compute( me, thread, Rs, F, Potential, Virial )
    type(tData), intent(inout) :: me
    integer,     intent(in)    :: thread
    real(rb),    intent(in)    :: Rs(3,me%natoms)
    real(rb),    intent(inout) :: F(3,me%natoms), Potential, Virial

    integer  :: i, j, k, m, n, icell, jcell, npairs, itype, nlocal, ntotal, first, last
    real(rb) :: xRc2, Rc2, r2, invR2, Eij, Wij, icharge
    logical  :: include(0:me%maxpairs)
    integer  :: atom(me%maxpairs), index(me%natoms)
    real(rb) :: Ri(3), Rij(3), Fi(3), Fij(3)
    integer,  allocatable :: xlist(:)
    type(tModel), pointer :: model

    xRc2 = me%xRcSq*me%invL2
    Rc2 = me%RcSq*me%invL2

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
          icharge = me%charge(i,me%layer)
          Ri = Rs(:,i)
          Fi = zero
          xlist = index(me%excluded%item(me%excluded%first(i):me%excluded%last(i)))
          include(xlist) = .false.
          associate (pairModel => me%model(:,itype,me%layer))
            do m = k + 1, ntotal
              if (include(m)) then
                j = atom(m)
                call c_f_pointer( pairModel(me%atomType(j)), model )
                if (associated(model)) then
                  Rij = Ri - Rs(:,j)
                  Rij = Rij - anint(Rij)
                  r2 = sum(Rij*Rij)
                  if (r2 < xRc2) then
                    npairs = npairs + 1
                    neighbor%item(npairs) = j
                    if (r2 < Rc2) then
                      invR2 = me%invL2/r2
                      call compute_pair()
                      Potential = Potential + Eij
                      Virial = Virial + Wij
                      Fij = Wij*invR2*Rij
                      Fi = Fi + Fij
                      F(:,j) = F(:,j) - Fij
                    end if
                  end if
               end if
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

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      include "compute_pair.f90"
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine find_pairs_and_compute

!===================================================================================================

  subroutine compute_pairs( me, thread, Rs, F, Potential, Virial )
    type(tData), intent(in)    :: me
    integer,     intent(in)    :: thread
    real(rb),    intent(in)    :: Rs(3,me%natoms)
    real(rb),    intent(inout) :: F(3,me%natoms), Potential, Virial

    integer  :: i, j, k, m, itype, firstAtom, lastAtom
    real(rb) :: Rc2, r2, invR2, Eij, Wij, icharge
    real(rb) :: Rij(3), Ri(3), Fi(3), Fij(3)
    type(tModel), pointer :: model

    Rc2 = me%RcSq*me%invL2
    firstAtom = me%cellAtom%first(me%threadCell%first(thread))
    lastAtom = me%cellAtom%last(me%threadCell%last(thread))
    associate (neighbor => me%neighbor(thread))
      do m = firstAtom, lastAtom
        i = me%cellAtom%item(m)
        itype = me%atomType(i)
        Ri = Rs(:,i)
        Fi = zero
        icharge = me%charge(i,me%layer)
        associate (pairModel => me%model(:,itype,me%layer))
          do k = neighbor%first(i), neighbor%last(i)
            j = neighbor%item(k)
            Rij = Ri - Rs(:,j)
            Rij = Rij - anint(Rij)
            r2 = sum(Rij*Rij)
            if (r2 < Rc2) then
              invR2 = me%invL2/r2
              call c_f_pointer( pairModel(me%atomType(j)), model )
              call compute_pair()
              Potential = Potential + Eij
              Virial = Virial + Wij
              Fij = Wij*invR2*Rij
              Fi = Fi + Fij
              F(:,j) = F(:,j) - Fij
            end if
          end do
        end associate
        F(:,i) = F(:,i) + Fi
      end do
    end associate

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      include "compute_pair.f90"
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine compute_pairs

!===================================================================================================

  subroutine compute_group_energy( me, thread, na, atom, energy )
    type(tData), intent(in)  :: me
    integer,     intent(in)  :: thread, na, atom(na)
    real(rb),    intent(out) :: energy(me%nlayers)

    integer  :: i, j, k, m, layer, itype, firstCell, lastCell
    real(rb) :: r2, invR2, Eij, Wij, icharge
    real(rb) :: Rij(3), Ri(3)
    type(tModel), pointer :: model

    firstCell = me%threadCell%first(thread)
    lastCell  = me%threadCell%last(thread)

    energy = zero
    associate (neighbor => me%neighbor(thread))
      do m = 1, size(atom)
        i = atom(m)
        if ((me%atomCell(i) >= firstCell).and.(me%atomCell(i) <= lastCell)) then
          itype = me%atomType(i)
          Ri = me%R(:,i)
          icharge = me%charge(i,me%layer)
          do k = neighbor%first(i), neighbor%last(i)
            j = neighbor%item(k)
            Rij = Ri - me%R(:,j)
            Rij = Rij - me%Lbox*anint(Rij*me%invL)
            r2 = sum(Rij*Rij)
            if (r2 < me%RcSq) then
              invR2 = one/r2
              do layer = 1, me%nlayers
                call c_f_pointer( me%model(me%atomType(j),itype,layer), model )
                call compute_pair()
                energy(layer) = energy(layer) + Eij
              end do
            end if
          end do
        end if
      end do
    end associate

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      include "compute_pair.f90"
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine compute_group_energy

!===================================================================================================

  subroutine EmDee_Rotational_Energies( md, Kr ) bind(C,name="EmDee_Rotational_Energies")
    type(tEmDee), value   :: md
    real(rb), intent(out) :: Kr(3)

    integer :: i
    type(tData), pointer :: me

    call c_f_pointer( md%data, me )

    Kr = zero
    do i = 1, me%nbodies
      Kr = Kr + me%body(i)%MoI*me%body(i)%omega**2
    end do
    Kr = half*Kr

  end subroutine EmDee_Rotational_Energies

!===================================================================================================

end module EmDeeCode
