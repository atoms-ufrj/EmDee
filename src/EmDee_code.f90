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

module EmDee_code

use omp_lib
use c_binding
use lists
use math
use models
use bonded_structs
use ArBee

implicit none

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

  integer(ib) :: builds         ! Number of neighbor-list builds
  real(rb)    :: time           ! Total time taken in force calculations
  real(rb)    :: Potential      ! Total potential energy of the system
  real(rb)    :: Kinetic        ! Total kinetic energy of the system
  real(rb)    :: Rotational     ! Rotational kinetic energy of the system
  real(rb)    :: Virial         ! Total internal virial of the system

  real(rb)    :: Lbox           ! Length of the simulation box
  type(c_ptr) :: coords         ! Pointer to the coordinates of all atoms
  type(c_ptr) :: momenta        ! Pointer to the momenta of all atoms
  type(c_ptr) :: forces         ! Pointer to the resultant forces on all atoms
  type(c_ptr) :: charge         ! Pointer to the electric charges of all atoms

  real(rb)    :: Rc             ! Cut-off distance
  real(rb)    :: RcSq           ! Cut-off distance squared
  real(rb)    :: xRc            ! Extended cutoff distance (including skin)
  real(rb)    :: xRcSq          ! Extended cutoff distance squared
  real(rb)    :: skinSq         ! Square of the neighbor list skin width
  real(rb)    :: invL           ! Inverse length of the simulation box
  real(rb)    :: invL2          ! Squared inverse length of the simulation box
  real(rb)    :: totalMass      ! Sum of the masses of all atoms
  integer(ib) :: coulomb        ! Flag for coulombic interactions
  real(rb)    :: eshift         ! Potential shifting factor for Coulombic interactions
  real(rb)    :: fshift         ! Force shifting factor for Coulombic interactions

  integer(ib) :: mcells         ! Number of cells at each dimension
  integer(ib) :: ncells         ! Total number of cells
  integer(ib) :: maxcells       ! Maximum number of cells
  integer(ib) :: maxatoms       ! Maximum number of atoms in a cell
  integer(ib) :: maxpairs       ! Maximum number of pairs formed by all atoms of a cell
  type(c_ptr) :: cell           ! Array containing all neighbor cells of each cell

  integer(ib) :: natoms         ! Number of atoms in the system
  type(c_ptr) :: atomType       ! Pointer to the type indexes of all atoms
  type(c_ptr) :: atomMass       ! Pointer to the masses of all atoms
  type(c_ptr) :: R0             ! Position of each atom at the latest neighbor list building

  integer(ib) :: ntypes         ! Number of atom types
  type(c_ptr) :: pairModel      ! Model of each type of atom pair

  type(c_ptr) :: bonds          ! List of bonds
  type(c_ptr) :: angles         ! List of angles
  type(c_ptr) :: dihedrals      ! List of dihedrals

  integer(ib) :: nbodies        ! Number of rigid bodies
  integer(ib) :: maxbodies      ! Maximum number of rigid bodies
  type(c_ptr) :: body           ! Pointer to the rigid bodies present in the system

  integer(ib) :: nfree          ! Number of independent atoms
  type(c_ptr) :: freeAtom       ! Pointer to the list of independent atoms

  integer(ib) :: nthreads       ! Number of parallel openmp threads
  integer(ib) :: threadAtoms    ! Number of atoms per parallel thread
  integer(ib) :: threadBodies   ! Number of rigid bodies per parallel thread
  type(c_ptr) :: cellAtom       ! List of atoms belonging to each cell
  type(c_ptr) :: threadCell     ! List of cells to be dealt with in each parallel thread
  type(c_ptr) :: neighbor       ! Pointer to neighbor lists
  type(c_ptr) :: excluded       ! List of pairs excluded from the neighbor lists
  type(c_ptr) :: random         ! Pointer for random number generators

end type tEmDee

private :: maximum_approach_sq, distribute_atoms, find_pairs_and_compute, &
           compute_pairs, compute_bonds, compute_angles, compute_dihedrals

contains

!===================================================================================================
!                                L I B R A R Y   P R O C E D U R E S
!===================================================================================================

  function EmDee_system( threads, rc, skin, N, types, masses, seed ) result( me ) &
                                                                     bind(C,name="EmDee_system")
    integer(ib), value :: threads, N, seed
    real(rb),    value :: rc, skin
    type(c_ptr), value :: types, masses
    type(tEmDee)       :: me

    integer(ib) :: i

    integer(ib),     pointer :: type_ptr(:)
    real(rb),        pointer :: mass_ptr(:)
    type(model_ptr), pointer :: pairModel(:,:)
    type(tList),     pointer :: cellAtom, threadCell, neighbor(:), excluded
    type(kiss),      pointer :: random(:)

    ! Set up fixed entities:
    me%nthreads = threads
    me%Rc = rc
    me%RcSq = rc*rc
    me%xRc = rc + skin
    me%xRcSq = me%xRc**2
    me%skinSq = skin*skin
    me%natoms = N
    me%fshift = one/me%RcSq
    me%eshift = -two/me%Rc

    if (c_associated(types)) then
      call c_f_pointer( types, type_ptr, [N] )
      if (minval(type_ptr) /= 1) stop "ERROR: wrong specification of atom types."
      me%ntypes = maxval(type_ptr)
      me%atomType = malloc_int( N, array = type_ptr )
    else
      me%ntypes = 1
      me%atomType = malloc_int( N, value = 1 )
    end if

    if (c_associated(masses)) then
      call c_f_pointer( masses, mass_ptr, [me%ntypes] )
      me%atomMass = malloc_real( N, array = mass_ptr(type_ptr) )
      me%totalMass = sum(mass_ptr(type_ptr))
    else
      me%atomMass = malloc_real( N, value = 1.0_rb )
      me%totalMass = real(N,rb)
    end if

    ! Initialize counters and other mutable entities:
    me%mcells = 0
    me%ncells = 0
    me%maxcells = 0
    me%builds = 0
    me%time = zero
    me%coulomb = 0
    me%Lbox = zero
    me%invL = huge(one)
    me%invL2 = huge(one)
    me%Potential = zero
    me%Kinetic = zero
    me%Rotational = zero
    me%coords = malloc_real( 3*N, value = zero )
    me%momenta = malloc_real( 3*N, value = zero )
    me%forces = malloc_real( 3*N, value = zero )
    me%charge = malloc_real( N, value = zero )
    me%R0 = malloc_real( 3*N, value = zero )
    me%cell = malloc_int( 0 )
    me%bonds = c_null_ptr
    me%angles = c_null_ptr
    me%dihedrals = c_null_ptr

    ! Allocate variables associated to rigid bodies:
    me%nbodies = 0
    me%maxbodies = 0
    me%body = c_null_ptr
    me%nfree = N
    me%freeAtom = malloc_int( N, array = [(i,i=1,N)] )
    me%threadAtoms = (N + threads - 1)/threads
    me%threadBodies = 0

    ! Allocate memory for list of atoms per cell:
    allocate( cellAtom )
    me%cellAtom = c_loc( cellAtom )
    call cellAtom % allocate( N, 0 )

    ! Allocate memory for lists of cells per parallel thread:
    allocate( threadCell )
    me%threadCell = c_loc( threadCell )

    ! Allocate memory for neighbor lists:
    allocate( neighbor(threads) )
    me%neighbor = c_loc(neighbor)
    call neighbor % allocate( extra, N )

    ! Allocate memory for the list of pairs excluded from the neighbor lists:
    allocate( excluded )
    me%excluded = c_loc( excluded )
    call excluded % allocate( extra, N )

    ! Allocate memory for pair models:
    allocate( pairModel(me%ntypes,me%ntypes) )
    me%pairModel = c_loc(pairModel(1,1))

    ! Initialize random number generators:
    allocate( random(threads) )
    me%random = c_loc( random )
    call random(1) % setup( seed )
    do i = 2, threads
      call random(i) % setup( random(i-1) % i32() )
    end do

  end function EmDee_system

!---------------------------------------------------------------------------------------------------
! TODO: copy charges from external array

  subroutine EmDee_set_charges( md, charges ) bind(C,name="EmDee_set_charges")
    type(c_ptr), value :: md, charges

    integer :: i, j

    type(tEmDee),      pointer :: me
    real(rb),          pointer :: chargesPtr(:)
    type(model_ptr),   pointer :: pairModel(:,:)
    type(EmDee_Model), pointer :: model

    call c_f_pointer( md, me )
    if (me%coulomb == 0) then
      call c_f_pointer( me%charge, chargesPtr, [me%natoms] )
      deallocate( chargesPtr )
      me%coulomb = 1
    end if
    me%charge = charges
    call c_f_pointer( me%pairModel, pairModel, [me%ntypes,me%ntypes] )
    do i = 1, me%ntypes
      do j = 1, me%ntypes
        model => pairModel(i,j)%model
        if (associated(model)) model%id = mCOULOMB + mod(model%id,mCOULOMB)
      end do
    end do

  end subroutine EmDee_set_charges

!---------------------------------------------------------------------------------------------------

  subroutine EmDee_set_pair_type( md, itype, jtype, model ) bind(C,name="EmDee_set_pair_type")
    type(c_ptr), value :: md
    integer(ib), value :: itype, jtype
    type(c_ptr), value :: model

    integer(ib) :: k
    logical :: keep

    type(tEmDee),      pointer :: me
    type(model_ptr),   pointer :: pairModel(:,:)
    type(EmDee_Model), pointer :: modelPtr, crossModel

    call c_f_pointer( md, me )
    call c_f_pointer( me%pairModel, pairModel, [me%ntypes,me%ntypes] )
    call c_f_pointer( model, modelPtr )
    if (itype == jtype) then
      call associate_model( itype, itype, modelPtr )
      do k = 1, me%ntypes
        if (k /= itype) then
          keep = associated(pairModel(itype,k)%model)
          if (keep) keep = pairModel(itype,k)%model%external == 1
          if (.not.keep) then
            crossModel => cross_pair( modelPtr, pairModel(k,k)%model )
            call replace( itype, k, crossModel )
          end if
        end if
      end do
    else
      call replace( itype, jtype, modelPtr )
    end if

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine associate_model( i, j, model )
        integer(ib),                intent(in)    :: i, j
        type(EmDee_Model), pointer, intent(inout) :: model
        if (associated(model)) then
          if (me%coulomb == 1) model%id = mCOULOMB + mod(model%id,mCOULOMB)
          pairModel(i,j)%model => model
        else
          nullify( pairModel(i,j)%model )
        end if
      end subroutine associate_model
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine replace( i, j, model )
        integer(ib),                intent(in)    :: i, j
        type(EmDee_Model), pointer, intent(inout) :: model
        type(EmDee_Model), pointer :: ijmodel
        ijmodel => pairModel(i,j)%model
        if (associated(ijmodel)) then
          if (ijmodel%external == 0) deallocate( ijmodel )
        end if
        call associate_model( i, j, model )
        call associate_model( j, i, model )
      end subroutine replace
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine EmDee_set_pair_type

!---------------------------------------------------------------------------------------------------

  subroutine EmDee_ignore_pair( md, i, j ) bind(C,name="EmDee_ignore_pair")
    type(c_ptr), value :: md
    integer(ib), value :: i, j

    integer(ib) :: n
    type(tEmDee), pointer :: me
    type(tList),  pointer :: excluded

    call c_f_pointer( md, me )
    if ((i > 0).and.(i <= me%natoms).and.(j > 0).and.(j <= me%natoms).and.(i /= j)) then
      call c_f_pointer( me%excluded, excluded )
      n = excluded%count
      if (n == excluded%nitems) call excluded % resize( n + extra )
      call add_item( i, j )
      call add_item( j, i )
      excluded%count = n
    end if

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine add_item( i, j )
        integer, intent(in) :: i, j
        integer :: start, end
        start = excluded%first(i)
        end = excluded%last(i)
        if ((end < start).or.(j > excluded%item(end))) then
          excluded%item(end+2:n+1) = excluded%item(end+1:n)
          excluded%item(end+1) = j
        else
          do while (j > excluded%item(start))
            start = start + 1
          end do
          if (j == excluded%item(start)) return
          excluded%item(start+1:n+1) = excluded%item(start:n)
          excluded%item(start) = j
        end if
        excluded%first(i+1:) = excluded%first(i+1:) + 1
        excluded%last(i:) = excluded%last(i:) + 1
        n = n + 1
      end subroutine add_item
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine EmDee_ignore_pair

!---------------------------------------------------------------------------------------------------

  subroutine EmDee_add_bond( md, i, j, model ) bind(C,name="EmDee_add_bond")
    type(c_ptr), value :: md
    integer(ib), value :: i, j
    type(c_ptr), value :: model

    type(tEmDee), pointer :: me

    call c_f_pointer( md, me )
    call add_bonded_struc( me%bonds, i, j, 0, 0, model )
    call EmDee_ignore_pair( md, i, j )

  end subroutine EmDee_add_bond

!---------------------------------------------------------------------------------------------------

  subroutine EmDee_add_angle( md, i, j, k, model ) bind(C,name="EmDee_add_angle")
    type(c_ptr), value :: md
    integer(ib), value :: i, j, k
    type(c_ptr), value :: model

    type(tEmDee), pointer :: me

    call c_f_pointer( md, me )
    call add_bonded_struc( me%angles, i, j, k, 0, model )
    call EmDee_ignore_pair( md, i, j )
    call EmDee_ignore_pair( md, i, k )
    call EmDee_ignore_pair( md, j, k )

  end subroutine EmDee_add_angle

!---------------------------------------------------------------------------------------------------

  subroutine EmDee_add_dihedral( md, i, j, k, l, model ) bind(C,name="EmDee_add_dihedral")
    type(c_ptr), value :: md
    integer(ib), value :: i, j, k, l
    type(c_ptr), value :: model

    type(tEmDee), pointer :: me

    call c_f_pointer( md, me )
    call add_bonded_struc( me%dihedrals, i, j, k, l, model )
    call EmDee_ignore_pair( md, i, j )
    call EmDee_ignore_pair( md, i, k )
    call EmDee_ignore_pair( md, i, l )
    call EmDee_ignore_pair( md, j, k )
    call EmDee_ignore_pair( md, j, l )
    call EmDee_ignore_pair( md, k, l )

  end subroutine EmDee_add_dihedral

!---------------------------------------------------------------------------------------------------

  subroutine EmDee_add_rigid_body( md, indexes, NP ) bind(C,name="EmDee_add_rigid_body")
    type(c_ptr), value :: md, indexes
    integer(ib), value :: NP

    integer  :: i, j

    logical, allocatable :: isolated(:)

    type(tEmDee),    pointer :: me
    integer(ib),     pointer :: atom(:), freeAtom(:)
    real(rb),        pointer :: R(:,:), mass(:)
    type(rigidBody), pointer :: body(:)

    call c_f_pointer( indexes, atom, [NP] )
    call c_f_pointer( md, me )
    call c_f_pointer( me%freeAtom, freeAtom, [me%nfree] )
    call c_f_pointer( me%atomMass, mass, [me%natoms] )
    call c_f_pointer( me%coords, R, [3,me%natoms] )
    if (me%nbodies == me%maxbodies) call realloc_rigid_body_list( me%body, me%maxbodies )
    call c_f_pointer( me%body, body, [me%nbodies+1] )

    allocate( isolated(me%nfree) )
    isolated = .true.
    do i = 1, NP
      where (freeAtom == atom(i)) isolated = .false.
    end do
    if (count(.not.isolated) /= NP) stop "Error adding rigid body: only free atoms are allowed."
    me%nfree = me%nfree - NP
    me%threadAtoms = (me%nfree + me%nthreads - 1)/me%nthreads
    freeAtom(1:me%nfree) = pack(freeAtom,isolated)
    do i = 1, NP-1
      do j = i+1, NP
        call EmDee_ignore_pair( md, atom(i), atom(j) )
      end do
    end do
    me%nbodies = me%nbodies + 1
    me%threadBodies = (me%nbodies + me%nthreads - 1)/me%nthreads
    call body(me%nbodies) % setup( atom, mass(atom) )
    if (me%Lbox /= zero) then
!TODO:    call body(me%nbodies) % update( )
    end if
!print*, body(me%nbodies) % q


  end subroutine EmDee_add_rigid_body

!---------------------------------------------------------------------------------------------------

  subroutine EmDee_set_coordinates( md, coords, L ) bind(C,name="EmDee_set_coordinates")
    type(c_ptr), value :: md, coords, L

    type(tEmDee),    pointer :: me
    integer(ib),     pointer :: free(:)
    real(rb),        pointer :: Lbox, Rext(:,:), Rsys(:,:)
    type(rigidBody), pointer :: body(:)

    call c_f_pointer( md, me )

    if (c_associated(L)) then
      call c_f_pointer( L, Lbox )
      me%Lbox = Lbox
      me%invL = one/Lbox
      me%invL2 = me%invL**2
    else if (me%Lbox == zero) then
      stop "ERROR: box side length has not been defined."
    end if

    if (c_associated(coords)) then
      call c_f_pointer( coords, Rext, [3,me%natoms] )
      call c_f_pointer( coords, Rsys, [3,me%natoms] )
      call c_f_pointer( me%freeAtom, free, [me%nfree] )
      call c_f_pointer( me%body, body, [me%nbodies] )
      !$omp parallel num_threads(me%nthreads)
      call assign_coordinates( omp_get_thread_num() + 1 )
      !$omp end parallel
    end if

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine assign_coordinates( thread )
        integer(ib), intent(in) :: thread
        integer(ib) :: i, j
        real(rb), allocatable :: R(:,:)
        do j = (thread - 1)*me%threadAtoms + 1, min(thread*me%threadAtoms, me%nfree)
          i = free(j)
          Rsys(:,i) = Rext(:,i)
        end do
        do i = (thread - 1)*me%threadBodies + 1, min(thread*me%threadBodies, me%nbodies)
          associate(b => body(i))
            R = Rext(:,b%index)
            forall (j=2:b%NP) R(:,j) = R(:,j) - me%Lbox*anint((R(:,j) - R(:,1))*me%invL)
            call b % update( R )
            Rsys(:,b%index) = R
          end associate
        end do
      end subroutine assign_coordinates
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine EmDee_set_coordinates

!---------------------------------------------------------------------------------------------------

  subroutine EmDee_generate_momenta( md, kT, adjust ) bind(C,name="EmDee_generate_momenta")
    type(c_ptr), value :: md
    real(rb),    value :: kT
    integer(ib), value :: adjust

    real(rb) :: KEt, KEr, pcm(3)

    real(rb),        pointer :: mass(:), P(:,:)
    type(tEmDee),    pointer :: me
    type(kiss),      pointer :: random(:)
    integer(ib),     pointer :: free(:)
    type(rigidBody), pointer :: body(:)

    call c_f_pointer( md, me )
    call c_f_pointer( me%atomMass, mass, [me%natoms] )
    call c_f_pointer( me%momenta, P, [3,me%natoms] )
    call c_f_pointer( me%random, random, [me%nthreads] )
    call c_f_pointer( me%body, body, [me%nbodies] )
    call c_f_pointer( me%freeAtom, free, [me%nbodies] )

    KEt = zero
    KEr = zero
    !$omp parallel num_threads(me%nthreads) reduction(+:KEt,KEr)
    call assign_momenta( omp_get_thread_num() + 1, KEt, KEr )
    !$omp end parallel
    me%Rotational = half*KEr
    me%Kinetic = half*(KEt + KEr)

    if (adjust == 1) then
      ! TODO: adjust kinetic energy and zero total momentum
    end if

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine assign_momenta( thread, KEt, KEr )
        integer(ib), intent(in)    :: thread
        real(rb),    intent(inout) :: KEt, KEr
        integer(ib) :: i, j
        real(rb) :: omega(3)
        associate (rng => random(thread))
          do j = (thread - 1)*me%threadAtoms + 1, min(thread*me%threadAtoms, me%nfree)
            i = free(j)
            P(:,i) = sqrt(mass(i)*kT)*[rng%normal(), rng%normal(), rng%normal()]
            KEt = KEt + sum(P(:,i)**2)/mass(i)
          end do
          do i = (thread - 1)*me%threadBodies + 1, min(thread*me%threadBodies, me%nbodies)
            associate (b => body(i))
              pcm = sqrt(b%mass*kT)*[rng%normal(), rng%normal(), rng%normal()]
              omega = sqrt(b%invMoI*kT)*[rng%normal(), rng%normal(), rng%normal()]
              call b % set_momenta( pcm, omega )
              P(:,b%index) = b%P
              KEt = KEt + b%invMass*sum(pcm*pcm)
              KEr = KEr + sum(b%MoI*omega*omega)
            end associate
          end do
        end associate
      end subroutine assign_momenta
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine EmDee_generate_momenta

!---------------------------------------------------------------------------------------------------

  subroutine EmDee_compute( md ) bind(C,name="EmDee_compute")
    type(c_ptr), value :: md

    integer(ib) :: M
    real(rb)    :: Potential, Virial
    logical     :: buildList

    type(tEmDee), pointer :: me
    real(rb),     pointer :: R(:,:), F(:,:), R0(:,:)

    real(rb), allocatable :: Rs(:,:), Fs(:,:)

    call c_f_pointer( md, me )
    me%time = me%time - omp_get_wtime()

    call c_f_pointer( me%forces, F, [3,me%natoms] )
    call c_f_pointer( me%coords, R, [3,me%natoms] )
    call c_f_pointer( me%R0, R0, [3,me%natoms] )

    allocate( Rs(3,me%natoms), Fs(3,me%natoms) )
    Rs = me%invL*R
    Fs = zero
    Potential = zero
    Virial = zero

    buildList = maximum_approach_sq( me%natoms, R - R0 ) > me%skinSq
    if (buildList) then
      M = floor(ndiv*me%Lbox/me%xRc)
      if (M < 5) stop "ERROR: simulation box is too small."
      call distribute_atoms( me, M, Rs )
      R0 = R
      me%builds = me%builds + 1
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
      if (c_associated(me%bonds)) call compute_bonds( me, thread, Rs, Fs, Potential, Virial )
      if (c_associated(me%angles)) call compute_angles( me, thread, Rs, Fs, Potential, Virial )
      if (c_associated(me%dihedrals)) call compute_dihedrals( me, thread, Rs, Fs, Potential, Virial )
    end block
    !$omp end parallel

    F = me%Lbox*Fs
    me%Potential = Potential
    me%Virial = third*Virial

    me%time = me%time + omp_get_wtime()

  end subroutine EmDee_compute

!===================================================================================================
!                              A U X I L I A R Y   P R O C E D U R E S
!===================================================================================================

  real(rb) function maximum_approach_sq( N, delta )
    integer(ib), intent(in) :: N
    real(rb),    intent(in) :: delta(3,N)
 
    integer(ib) :: i
    real(rb)    :: maximum, next, deltaSq

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

!---------------------------------------------------------------------------------------------------

  subroutine distribute_atoms( me, M, Rs )
    type(tEmDee), intent(inout) :: me
    integer(ib),  intent(in)    :: M
    real(rb),     intent(in)    :: Rs(3,me%natoms)

    integer(ib) :: MM, cells_per_thread, maxNatoms
    logical     :: make_cells
    integer(ib) :: atomCell(me%natoms), threadNatoms(me%nthreads), next(me%natoms)
    integer(ib), allocatable :: natoms(:)

    type(tCell), pointer :: cell(:)
    type(tList), pointer :: threadCell, cellAtom

    call c_f_pointer( me%cell, cell, [me%maxcells] )
    call c_f_pointer( me%threadCell, threadCell )
    call c_f_pointer( me%cellAtom, cellAtom )

    MM = M*M
    make_cells = M /= me%mcells
    if (make_cells) then
      me%mcells = M
      me%ncells = M*MM
      if (me%ncells > me%maxcells) then
        deallocate( cell, cellAtom%first, cellAtom%last )
        allocate( cell(me%ncells), cellAtom%first(me%ncells), cellAtom%last(me%ncells) )
        call threadCell % allocate( 0, me%nthreads )
        me%maxcells = me%ncells
        me%cell = c_loc(cell(1))
      end if
      cells_per_thread = (me%ncells + me%nthreads - 1)/me%nthreads
    end if

    allocate( natoms(me%ncells) )

    !$omp parallel num_threads(me%nthreads) reduction(max:maxNatoms)
    block
      integer(ib) :: thread, i, k, icell, ix, iy, iz, first, last, atoms_per_thread
      integer(ib) :: icoord(3)
      integer(ib), allocatable :: head(:)

      thread = omp_get_thread_num() + 1

      if (make_cells) then
        first = (thread - 1)*cells_per_thread + 1
        last = min( thread*cells_per_thread, me%ncells )
        do icell = first, last
          k = icell - 1
          iz = k/MM
          iy = (k - iz*MM)/M
          ix = k - (iy*M + iz*MM)
          cell(icell)%neighbor = 1 + pbc(ix+nb(1,:)) + pbc(iy+nb(2,:))*M + pbc(iz+nb(3,:))*MM
        end do
        threadCell%first(thread) = first
        threadCell%last(thread) = last
      else
        first = threadCell%first(thread)
        last = threadCell%last(thread)
      end if

      atoms_per_thread = (me%natoms + me%nthreads - 1)/me%nthreads
      do i = (thread - 1)*atoms_per_thread + 1, min( thread*atoms_per_thread, me%natoms )
        icoord = int(M*(Rs(:,i) - floor(Rs(:,i))),ib)
        atomCell(i) = 1 + icoord(1) + M*icoord(2) + MM*icoord(3)
      end do
      !$omp barrier

      allocate( head(first:last) )
      head = 0
      natoms(first:last) = 0
      do i = 1, me%natoms
        icell = atomCell(i)
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
        cellAtom%first(icell) = k + 1
        i = head(icell)
        do while (i /= 0)
          k = k + 1
          cellAtom%item(k) = i
          i = next(i)
        end do
        cellAtom%last(icell) = k
        if (natoms(icell) > maxNatoms) maxNatoms = natoms(icell)
      end do

    end block
    !$omp end parallel
    me%maxatoms = maxNatoms
    me%maxpairs = (maxNatoms*((2*nbcells + 1)*maxNatoms - 1))/2

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      elemental integer(ib) function pbc( x )
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
  end subroutine distribute_atoms

!---------------------------------------------------------------------------------------------------

  subroutine compute_bonds( me, threadId, R, F, Potential, Virial )
    type(tEmDee), intent(inout) :: me
    integer,      intent(in)    :: threadId
    real(rb),     intent(in)    :: R(3,me%natoms)
    real(rb),     intent(inout) :: F(3,me%natoms), Potential, Virial

    integer(ib) :: i, j, m, nbonds
    real(rb)    :: d, E, mdEdr
    real(rb)    :: Rij(3), Fij(3)

    type(EmDee_Model), pointer :: model
    type(tStructData), pointer :: bonds

    call c_f_pointer( me%bonds, bonds )

    nbonds = (bonds%number + me%nthreads - 1)/me%nthreads
    do m = (threadId - 1)*nbonds + 1, min( bonds%number, threadId*nbonds )
      i = bonds%item(m)%i
      j = bonds%item(m)%j
      Rij = R(:,i) - R(:,j)
      Rij = Rij - anint(Rij)
      d = me%Lbox*sqrt(sum(Rij*Rij))
      model => bonds%item(m)%model
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

!---------------------------------------------------------------------------------------------------

  subroutine compute_angles( me, threadId, R, F, Potential, Virial )
    type(tEmDee), intent(inout) :: me
    integer,      intent(in)    :: threadId
    real(rb),     intent(in)    :: R(3,me%natoms)
    real(rb),     intent(inout) :: F(3,me%natoms), Potential, Virial

    integer(ib) :: i, j, k, m, nangles
    real(rb)    :: aa, bb, ab, axb, theta, Ea, Fa
    real(rb)    :: Rj(3), Fi(3), Fk(3), a(3), b(3)

    type(EmDee_Model), pointer :: model
    type(tStructData), pointer :: angles

    call c_f_pointer( me%angles, angles )

    nangles = (angles%number + me%nthreads - 1)/me%nthreads
    do m = (threadId - 1)*nangles + 1, min( angles%number, threadId*nangles )
      i = angles%item(m)%i
      j = angles%item(m)%j
      k = angles%item(m)%k
      model => angles%item(m)%model
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

!---------------------------------------------------------------------------------------------------

  subroutine compute_dihedrals( me, threadId, R, F, Potential, Virial )
    type(tEmDee), intent(inout) :: me
    integer,      intent(in)    :: threadId
    real(rb),     intent(in)    :: R(3,me%natoms)
    real(rb),     intent(inout) :: F(3,me%natoms), Potential, Virial

    integer(ib) :: m, ndihedrals, i, j
    real(rb)    :: Rc2, Ed, Fd, r2, invR2, Eij, Wij, icharge
    real(rb)    :: Rj(3), Rk(3), Fi(3), Fk(3), Fl(3), Fij(3)
    real(rb)    :: normRkj, normX, a, b, phi
    real(rb)    :: rij(3), rkj(3), rlk(3), x(3), y(3), z(3), u(3), v(3), w(3)

    integer(ib),       pointer :: atomType(:)
    real(rb),          pointer :: charge(:)
    type(tStruct),     pointer :: d
    type(tStructData), pointer :: dihedrals
    type(EmDee_Model), pointer :: model
    type(model_ptr),   pointer :: pairModel(:,:)

    call c_f_pointer( me%dihedrals, dihedrals )
    call c_f_pointer( me%charge, charge, [me%natoms] )
    call c_f_pointer( me%atomType, atomType, [me%natoms] )
    call c_f_pointer( me%pairModel, pairModel, [me%ntypes,me%ntypes] )

    Rc2 = me%RcSq*me%invL2
    ndihedrals = (dihedrals%number + me%nthreads - 1)/me%nthreads
    do m = (threadId - 1)*ndihedrals + 1, min( dihedrals%number, threadId*ndihedrals )
      d => dihedrals%item(m)
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
          model => pairModel(atomType(i),atomType(j))%model
          icharge = charge(i)
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

!---------------------------------------------------------------------------------------------------

  subroutine find_pairs_and_compute( me, thread, Rs, F, Potential, Virial )
    type(tEmDee), intent(inout) :: me
    integer,      intent(in)    :: thread
    real(rb),     intent(in)    :: Rs(3,me%natoms)
    real(rb),     intent(inout) :: F(3,me%natoms), Potential, Virial

    integer(ib) :: i, j, k, m, n, icell, jcell, npairs, itype, nlocal, ntotal, first, last
    real(rb)    :: xRc2, Rc2, r2, invR2, Eij, Wij, icharge
    logical     :: include(0:me%maxpairs)
    integer(ib) :: atom(me%maxpairs), index(me%natoms)
    real(rb)    :: Ri(3), Rij(3), Fi(3), Fij(3)

    integer(ib), allocatable :: xlist(:)

    type(tCell),       pointer :: cell(:)
    integer(ib),       pointer :: atomType(:)
    real(rb),          pointer :: charge(:)
    type(EmDee_Model), pointer :: model
    type(model_ptr),   pointer :: pairModel(:,:)
    type(tList),       pointer :: cellAtom, threadCell, neighborLists(:), neighbor, excluded

    call c_f_pointer( me%cell, cell, [me%ncells] )
    call c_f_pointer( me%atomType, atomType, [me%natoms] )
    call c_f_pointer( me%charge, charge, [me%natoms] )
    call c_f_pointer( me%pairModel, pairModel, [me%ntypes,me%ntypes] )
    call c_f_pointer( me%cellAtom, cellAtom )
    call c_f_pointer( me%threadCell, threadCell )
    call c_f_pointer( me%excluded, excluded )
    call c_f_pointer( me%neighbor, neighborLists, [me%nthreads] )
    neighbor => neighborLists(thread)

    xRc2 = me%xRcSq*me%invL2
    Rc2 = me%RcSq*me%invL2

    include = .true.
    index = 0
    npairs = 0
    do icell = threadCell%first(thread), threadCell%last(thread)

      if (neighbor%nitems < npairs + me%maxpairs) then
        call neighbor % resize( npairs + me%maxpairs + extra )
      end if

      first = cellAtom%first(icell)
      last = cellAtom%last(icell)
      nlocal = last - first + 1
      atom(1:nlocal) = cellAtom%item(first:last)

      ntotal = nlocal
      do m = 1, nbcells
        jcell = cell(icell)%neighbor(m)
        first = cellAtom%first(jcell)
        last = cellAtom%last(jcell)
        n = ntotal + 1
        ntotal = n + last - first
        atom(n:ntotal) = cellAtom%item(first:last)
      end do

      forall (m=1:ntotal) index(atom(m)) = m
      do k = 1, nlocal
        i = atom(k)
        neighbor%first(i) = npairs + 1
        itype = atomType(i)
        icharge = charge(i)
        Ri = Rs(:,i)
        Fi = zero
        xlist = index(excluded%item(excluded%first(i):excluded%last(i)))
        include(xlist) = .false.
        do m = k + 1, ntotal
          if (include(m)) then
            j = atom(m)
            model => pairModel(itype,atomType(j))%model
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
        F(:,i) = F(:,i) + Fi
        neighbor%last(i) = npairs
        include(xlist) = .true.
      end do      
      index(atom(1:ntotal)) = 0

    end do
    neighbor%count = npairs
    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      include "compute_pair.f90"
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine find_pairs_and_compute

!---------------------------------------------------------------------------------------------------

  subroutine compute_pairs( me, thread, Rs, F, Potential, Virial )
    type(tEmDee), intent(in)    :: me
    integer,      intent(in)    :: thread
    real(rb),     intent(in)    :: Rs(3,me%natoms)
    real(rb),     intent(inout) :: F(3,me%natoms), Potential, Virial

    integer  :: i, j, k, m, itype, firstAtom, lastAtom
    real(rb) :: Rc2, r2, invR2, Eij, Wij, icharge
    real(rb) :: Rij(3), Ri(3), Fi(3), Fij(3)

    integer(ib),       pointer :: atomType(:)
    real(rb),          pointer :: charge(:)
    type(EmDee_Model), pointer :: model
    type(model_ptr),   pointer :: pairModel(:,:)
    type(tList),       pointer :: cellAtom, threadCell, neighborLists(:), neighbor

    call c_f_pointer( me%atomType, atomType, [me%natoms] )
    call c_f_pointer( me%charge, charge, [me%natoms] )
    call c_f_pointer( me%pairModel, pairModel, [me%ntypes,me%ntypes] )
    call c_f_pointer( me%cellAtom, cellAtom )
    call c_f_pointer( me%threadCell, threadCell )
    call c_f_pointer( me%neighbor, neighborLists, [me%nthreads] )
    neighbor => neighborLists(thread)

    Rc2 = me%RcSq*me%invL2
    firstAtom = cellAtom%first(threadCell%first(thread))
    lastAtom = cellAtom%last(threadCell%last(thread))
    do m = firstAtom, lastAtom
      i = cellAtom%item(m)
      itype = atomType(i)
      Ri = Rs(:,i)
      Fi = zero
      icharge = charge(i)
      do k = neighbor%first(i), neighbor%last(i)
        j = neighbor%item(k)
        Rij = Ri - Rs(:,j)
        Rij = Rij - anint(Rij)
        r2 = sum(Rij*Rij)
        if (r2 < Rc2) then
          invR2 = me%invL2/r2
          model => pairModel(itype,atomType(j))%model
          call compute_pair()
          Potential = Potential + Eij
          Virial = Virial + Wij
          Fij = Wij*invR2*Rij
          Fi = Fi + Fij
          F(:,j) = F(:,j) - Fij
        end if
      end do
      F(:,i) = F(:,i) + Fi
    end do

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      include "compute_pair.f90"
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine compute_pairs

!---------------------------------------------------------------------------------------------------

end module EmDee_code
