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

! TODO: Check kspace discounts in rigid bodies (presence of kCoul)
! TODO: Discount long-range electrostatic terms of dihedrals
! TODO: Replace type(c_ptr) arguments by array arguments when possible
! TODO: Create indexing for having sequential body particles and free particles in arrays

module EmDeeCode

#define PRINT .false.

use EmDeeData
use neighbor_lists

implicit none

private

character(11), parameter :: VERSION = "10 Oct 2018"

type, bind(C), public :: tOpts
  logical(lb) :: Translate            ! Flag to activate/deactivate translations
  logical(lb) :: Rotate               ! Flag to activate/deactivate rotations
  integer(ib) :: RotationMode         ! Algorithm used for free rotation of rigid bodies
  logical(lb) :: AutoBodyUpdate       ! Flag to activate/deactivate automatic rigid body update
  logical(lb) :: Compute              ! Flag to activate/deactivate energy computations
end type tOpts

type, bind(C), public :: tTime
  real(rb) :: Pair
  real(rb) :: Motion
  real(rb) :: Neighbor
  real(rb) :: Total
end type tTime

type, bind(C), public :: tEmDee
  integer(ib)    :: Builds             ! Number of neighbor list builds
  type(tTime)    :: Time
  type(tEnergy)  :: Energy             ! All potential energy terms
  type(tKinetic) :: Kinetic            ! All kinetic energy terms
  type(tVirial)  :: Virial
  integer(ib)    :: DoF                ! Total number of degrees of freedom
  integer(ib)    :: RotDoF             ! Number of rotational degrees of freedom
  type(c_ptr)    :: Data               ! Pointer to system data
  type(tOpts)    :: Options            ! List of options to change EmDee's behavior
end type tEmDee

contains

!===================================================================================================
!                                L I B R A R Y   P R O C E D U R E S
!===================================================================================================

  function EmDee_system( threads, layers, rc, skin, N, types, masses, bodies ) &
    bind(C,name="EmDee_system")
    integer(ib), value :: threads, layers, N
    real(rb),    value :: rc, skin
    type(c_ptr), value :: types, masses, bodies
    type(tEmDee)       :: EmDee_system

    integer,     pointer :: ptype(:)
    real(rb),    pointer :: pmass(:)
    type(tData), pointer :: me
    type(pairContainer) :: noPair
    type(coulContainer) :: noCoul

    write(*,'("EmDee (version: ",A11,")")') VERSION

    ! Allocate data structure:
    allocate( me )

    ! Set up fixed entities:
    me%nthreads = threads
    me%nlayers = layers
    me%Rc = rc
    me%skin = skin
    me%RcSq = rc*rc
    me%xRc = rc + skin
    me%xRcSq = me%xRc**2
    me%InRc = me%Rc
    me%InRcSq = me%RcSq
    me%xInRcSq = me%xRcSq
    me%skinSq = skin*skin
    me%natoms = N
    me%threadAtoms = (N + threads - 1)/threads
    me%kspace_active = .false.

    ! Set up atom types:
    if (c_associated(types)) then
      call c_f_pointer( types, ptype, [N] )
      if (minval(ptype) /= 1) call error( "system setup", "wrong specification of atom types" )
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
    allocate( me%P(3,N), me%R0(3,N), source = zero )
    allocate( me%charge(N), source = zero )
    allocate( me%charged(N), source = .false. )
    allocate( me%cell(0), me%atomsInCell(0) )
    allocate( me%atomCell(N) )

    ! Allocate variables related to rigid bodies:
    call allocate_rigid_bodies( me, bodies )

    ! Allocate memory for list of atoms per cell:
    call me % cellAtom % allocate( N, 0 )

    ! Allocate memory for neighbor lists:
    allocate( me%neighbor(threads) )
    call me % neighbor % allocate( extra, N, middle = .true., value = .true. )

    ! Allocate memory for the list of pairs excluded from the neighbor lists:
    call me % excluded % allocate( extra, N )

    ! Allocate memory for variables regarding pair models:
    allocate( pair_none :: noPair%model )
    call noPair % model % setup()
    allocate( me%pair(me%ntypes,me%ntypes,me%nlayers), source = noPair )
    allocate( me%multilayer(me%ntypes,me%ntypes), source = .false. )
    allocate( me%overridable(me%ntypes,me%ntypes), source = .true. )
    allocate( me%interact(me%ntypes,me%ntypes), source = .false. )

    ! Allocate memory for variables regarding Coulomb models:
    allocate( coul_none :: noCoul%model )
    call noCoul % model % setup()
    me%multilayer_coulomb = .false.
    allocate( me%coul(me%nlayers), source = noCoul )

    ! Allocate other variables related to model layers:
    allocate( me%pairs_exist(me%nlayers), source = .false. )
    allocate( me%useInRc(layers), source = .false. )
    allocate( me%bonded(layers), source = .true. )
    allocate( me%forcesUpToDate(layers), source = .false. )
    allocate( me%layerF(3,N,layers), source = zero )
    me%layer = 1
    me%F => me%layerF(:,:,1)

    ! Set system as uninitialized:
    me%initialized = .false.

    ! Set up mutable entities:
    EmDee_system % builds = 0

    EmDee_system % Time % Pair = zero
    EmDee_system % Time % Motion = zero
    EmDee_system % Time % Neighbor = zero
    EmDee_system % Time % Total = zero

    EmDee_system % Energy % Potential = zero
    EmDee_system % Energy % Dispersion = zero
    EmDee_system % Energy % Coulomb = zero
    EmDee_system % Kinetic % Total = zero
    EmDee_system % Kinetic % TransPart = zero
    EmDee_system % Kinetic % Rotational = zero
    EmDee_system % Kinetic % RotPart = zero
    EmDee_system % Energy % ShadowPotential = zero
    EmDee_system % Kinetic % ShadowKinetic = zero
    EmDee_system % Kinetic % ShadowRotational = zero
    EmDee_system % Energy % UpToDate = .false.

    allocate( me%layerEnergy(layers), source = EmDee_system % Energy )
    allocate( me%layerVirial(layers), source = EmDee_system % Virial )

    EmDee_system % Virial % Total = zero
    EmDee_system % Virial % Body = zero
    EmDee_system % DoF = 3*(N - 1)
    EmDee_system % RotDoF = 0
    EmDee_system % data = c_loc(me)
    EmDee_system % Options % translate = .true.

    EmDee_system % Options % Rotate = .true.
    EmDee_system % Options % RotationMode = 0
    EmDee_system % Options % AutoBodyUpdate = .true.
    EmDee_system % Options % Compute = .true.

  end function EmDee_system

!===================================================================================================

  subroutine EmDee_share_phase_space( mdkeep, mdlose ) bind(C,name="EmDee_share_phase_space")
    type(tEmDee), value         :: mdkeep
    type(tEmDee), intent(inout) :: mdlose

    character(*), parameter :: task = "phase space sharing"

    type(tData), pointer :: keep, lose

    call c_f_pointer( mdkeep%data, keep )
    call c_f_pointer( mdlose%data, lose )

    if (.not.(keep%initialized.and.lose%initialized)) &
      call error( task, "EmDee system 1 has not been initialized" )

    if (keep%natoms /= lose%natoms) call error( task, "different numbers of atoms" )
    if (any(keep%atomType /= lose%atomType)) call error( task, "atom types do not match" )
    if (any(keep%mass /= lose%mass)) call error( task, "atom masses do not match" )
    if (any(keep%atomBody /= lose%atomBody)) call error( task, "rigid bodies do not match" )

    if (lose%initialized) deallocate( lose%R, lose%P, lose%body, lose%Lbox )
    lose%R => keep%R
    lose%P => keep%P
    lose%body => keep%body
    lose%Lbox => keep%Lbox

    mdlose%Kinetic%Total = mdkeep%Kinetic%Total
    mdlose%Kinetic%TransPart = mdkeep%Kinetic%TransPart
    mdlose%Kinetic%Rotational = mdkeep%Kinetic%Rotational
    mdlose%Kinetic%RotPart = mdkeep%Kinetic%RotPart

  end subroutine EmDee_share_phase_space

!===================================================================================================

  subroutine EmDee_layer_based_parameters( md, InternalRc, Apply, Bonded ) &
    bind(C,name="EmDee_layer_based_parameters")
    type(tEmDee), value      :: md
    real(rb),     value      :: InternalRc
    integer(ib),  intent(in) :: Apply(*)
    integer(ib),  intent(in) :: Bonded(*)

    character(*), parameter :: task = "layer-based parameter setting"

    integer :: layer
    type(tData), pointer :: me

    call c_f_pointer( md%data, me )

    if (me%initialized) call error(task, "system has already been initialized")

    me%useInRc = Apply(1:me%nlayers) /= 0
    if (any(me%useInRc).and.((InternalRc <= zero).or.(InternalRc > me%Rc))) then
      call error(task, "invalid internal cutoff specification")
    end if
    me%InRc = InternalRc
    me%InRcSq = InternalRc**2
    me%xInRcSq = (InternalRc + me%skin)**2

    ! Mark layers with bonded interactions:
    me%bonded = Bonded(1:me%nlayers) /= 0

    ! Reinitialize coulomb models:
    do layer = 1, me%nlayers
      call me % coul(layer) % model % cutoff_setup( merge(me%InRc, me%Rc, me%useInRc(layer)) )
    end do

  end subroutine EmDee_layer_based_parameters

!===================================================================================================

  subroutine EmDee_set_pair_model( md, itype, jtype, model, kCoul ) &
    bind(C,name="EmDee_set_pair_model")
    type(tEmDee), value :: md
    integer(ib),  value :: itype, jtype
    type(c_ptr),  value :: model
    real(rb),     value :: kCoul

    character(*), parameter :: task = "pair model setup"

    integer :: layer, ktype
    type(tData), pointer :: me
    type(modelContainer), pointer :: container

    call c_f_pointer( md%data, me )
    if (me%initialized) then
      call error( task, "cannot set model after coordinates have been defined" )
    end if

    if (.not.ranged( [itype,jtype], me%ntypes )) then
      call error( task, "provided type index is out of range" )
    end if

    if (.not.c_associated(model)) then
      call error( task, "a valid pair model must be provided" )
    end if
    call c_f_pointer( model, container )

    select type ( pair => container%model )
      class is (cPairModel)
        do layer = 1, me%nlayers
          call set_pair_type( me, itype, jtype, layer, container, kCoul )
        end do
      class default
        call error( task, "a valid pair model must be provided" )
    end select

    me%multilayer(itype,jtype) = .false.
    me%multilayer(jtype,itype) = .false.
    if (itype == jtype) then
      do ktype = 1, me%ntypes
        if ((ktype /= itype).and.me%overridable(itype,ktype)) then
          me%multilayer(itype,ktype) = me%multilayer(ktype,ktype)
          me%multilayer(ktype,itype) = me%multilayer(ktype,ktype)
        end if
      end do
    else
      me%overridable(itype,jtype) = .false.
      me%overridable(jtype,itype) = .false.
    end if

  end subroutine EmDee_set_pair_model

!===================================================================================================

  subroutine EmDee_set_pair_multimodel( md, itype, jtype, model, kCoul ) &
    bind(C,name="EmDee_set_pair_multimodel")
    type(tEmDee), value      :: md
    integer(ib),  value      :: itype, jtype
    type(c_ptr),  intent(in) :: model(*)
    real(rb),     intent(in) :: kCoul(*)

    character(*), parameter :: task = "pair multimodel setup"

    integer :: layer, ktype
    character(5) :: C
    type(tData), pointer :: me
    type(modelContainer), pointer :: container

    call c_f_pointer( md%data, me )
    if (me%initialized) call error( task, "cannot set model after coordinates have been defined" )

    if (.not.ranged( [itype,jtype], me%ntypes )) then
      call error( task, "provided type index is out of range" )
    end if

    do layer = 1, me%nlayers
      if (.not.c_associated(model(layer))) then
        write(C,'(I5)') me%nlayers
        call error( task, trim(adjustl(C))//" valid pair models must be provided" )
      end if
      call c_f_pointer( model(layer), container )
      select type ( pair => container%model )
        class is (cPairModel)
          call set_pair_type( me, itype, jtype, layer, container, kCoul(layer) )
        class default
          call error( task, "a valid pair model must be provided" )
      end select
    end do

    me%multilayer(itype,jtype) = .true.
    me%multilayer(jtype,itype) = .true.
    if (itype == jtype) then
      do ktype = 1, me%ntypes
        if ((ktype /= itype).and.me%overridable(itype,ktype)) then
          me%multilayer(itype,ktype) = .true.
          me%multilayer(ktype,itype) = .true.
        end if
      end do
    else
      me%overridable(itype,jtype) = .false.
      me%overridable(jtype,itype) = .false.
    end if

  end subroutine EmDee_set_pair_multimodel

!===================================================================================================

  subroutine EmDee_set_kspace_model( md, model ) bind(C,name="EmDee_set_kspace_model")
    type(tEmDee), value :: md
    type(c_ptr),  value :: model

    character(*), parameter :: task = "kspace model setup"

    type(tData), pointer :: me
    type(modelContainer), pointer :: container

    call c_f_pointer( md%data, me )
    if (me%initialized) then
      call error( task, "cannot set model after coordinates have been defined" )
    end if

    if (.not.c_associated(model)) then
      call error( task, "a valid kspace model must be provided" )
    end if
    call c_f_pointer( model, container )

    select type ( kspace => container%model )
      class is (cKspaceModel)
        if (me%kspace_active) deallocate( me%kspace )
        allocate( me%kspace, source = kspace )
        me%kspace_active = .true.
      class default
        call error( task, "a valid kspace model must be provided" )
    end select

  end subroutine EmDee_set_kspace_model

!===================================================================================================

  subroutine EmDee_set_coul_model( md, model ) bind(C,name="EmDee_set_coul_model")
    type(tEmDee), value :: md
    type(c_ptr),  value :: model

    character(*), parameter :: task = "coulomb model setup"

    integer  :: layer
    real(rb) :: layerRc
    type(tData), pointer :: me
    type(modelContainer), pointer :: container

    call c_f_pointer( md%data, me )
    if (me%initialized) then
      call error( task, "cannot set model after coordinates have been defined" )
    end if

    if (.not.c_associated(model)) then
      call error( task, "a valid coulomb model must be provided" )
    end if
    call c_f_pointer( model, container )

    select type ( coul => container%model )
      class is (cCoulModel)
        do layer = 1, me%nlayers
          layerRc = merge(me%InRc, me%Rc, me%useInRc(layer))
          me%coul(layer) = container
          call me % coul(layer) % model % cutoff_setup( layerRc )
          call me % coul(layer) % model % modifier_setup( layerRc )
        end do
      class default
        call error( task, "a valid coulomb model must be provided" )
    end select
    me%multilayer_coulomb = .false.

  end subroutine EmDee_set_coul_model

!===================================================================================================

  subroutine EmDee_set_coul_multimodel( md, model ) bind(C,name="EmDee_set_coul_multimodel")
    type(tEmDee), value :: md
    type(c_ptr),  intent(in) :: model(*)

    character(*), parameter :: task = "coulomb multimodel setup"

    integer  :: layer
    real(rb) :: layerRc
    character(5) :: C
    type(tData), pointer :: me
    type(modelContainer), pointer :: container

    call c_f_pointer( md%data, me )
    if (me%initialized) call error( task, "cannot set model after coordinates have been defined" )

    write(C,'(I5)') me%nlayers
    C = adjustl(C)
    do layer = 1, me%nlayers
      layerRc = merge(me%InRc, me%Rc, me%useInRc(layer))
      if (.not.c_associated(model(layer))) then
        call error( task, trim(C)//" valid coulomb models must be provided" )
      end if
      call c_f_pointer( model(layer), container )
      select type ( coul => container%model )
        class is (cCoulModel)
          me%coul(layer) = container
          call me % coul(layer) % model % cutoff_setup( layerRc )
          call me % coul(layer) % model % modifier_setup( layerRc )
        class default
          call error( task, "a valid coulomb model must be provided" )
      end select
    end do
    me%multilayer_coulomb = .true.

  end subroutine EmDee_set_coul_multimodel

!===================================================================================================

  subroutine EmDee_ignore_pair( md, i, j ) bind(C,name="EmDee_ignore_pair")
    type(tEmDee), value :: md
    integer(ib),  value :: i, j

    integer :: n
    type(tData), pointer :: me

    call c_f_pointer( md%data, me )
    if ((i /= j).and.ranged([i,j],me%natoms)) then
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
        if (end < start) then
          excluded%item(end+2:n+1) = excluded%item(end+1:n)
          excluded%item(end+1) = j
        elseif (j > excluded%item(end)) then ! Repetition avoids exception in DEBUB mode
          excluded%item(end+2:n+1) = excluded%item(end+1:n)
          excluded%item(end+1) = j
        else
          do while (j > excluded%item(start))
            start = start + 1
          end do
          if (j == excluded%item(start)) return
          excluded%item(start+1:n+1) = excluded%item(start:n)
          excluded%item(start) = j
          start = start + 1
        end if
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

    type(tData),          pointer :: me
    type(modelContainer), pointer :: container

    call c_f_pointer( md%data, me )

    if (.not.ranged([i,j],me%natoms)) call error( "add_bond", "atom index out of range" )
    if (.not.c_associated(model)) call error( "add_bond", "a valid model must be provided" )
    call c_f_pointer( model, container )

    select type (bmodel => container%model)
      class is (cBondModel)
        call me % bonds % add( i, j, 0, 0, bmodel )
        call EmDee_ignore_pair( md, i, j )
      class default
        call error( "add_bond", "the provided model must be a bond model" )
    end select

  end subroutine EmDee_add_bond

!===================================================================================================

  subroutine EmDee_add_angle( md, i, j, k, model ) bind(C,name="EmDee_add_angle")
    type(tEmDee), value :: md
    integer(ib),  value :: i, j, k
    type(c_ptr),  value :: model

    type(tData), pointer :: me
    type(modelContainer), pointer :: container

    call c_f_pointer( md%data, me )

    if (.not.ranged([i,j,k],me%natoms)) call error( "add_angle", "atom index out of range" )
    if (.not.c_associated(model)) call error( "add_angle", "a valid model must be provided" )
    call c_f_pointer( model, container )

    select type (amodel => container%model)
      class is (cAngleModel)
        call me % angles % add( i, j, k, 0, amodel )
        call EmDee_ignore_pair( md, i, j )
        call EmDee_ignore_pair( md, i, k )
        call EmDee_ignore_pair( md, j, k )
      class default
        call error( "add_angle", "the provided model must be an angle model" )
    end select

  end subroutine EmDee_add_angle

!===================================================================================================

  subroutine EmDee_add_dihedral( md, i, j, k, l, model ) bind(C,name="EmDee_add_dihedral")
    type(tEmDee), value :: md
    integer(ib),  value :: i, j, k, l
    type(c_ptr),  value :: model

    type(tData), pointer :: me
    type(modelContainer), pointer :: container

    call c_f_pointer( md%data, me )

    if (.not.ranged([i,j,k,l],me%natoms)) call error( "add_dihedral", "atom index out of range" )
    if (.not.c_associated(model)) call error( "add_dihedral", "a valid model must be provided" )
    call c_f_pointer( model, container )

    select type (dmodel => container%model)
      class is (cDihedralModel)
        call me % dihedrals % add( i, j, k, l, dmodel )
        call EmDee_ignore_pair( md, i, j )
        call EmDee_ignore_pair( md, i, k )
        call EmDee_ignore_pair( md, i, l )
        call EmDee_ignore_pair( md, j, k )
        call EmDee_ignore_pair( md, j, l )
        call EmDee_ignore_pair( md, k, l )
      class default
        call error( "add_dihedral", "the provided model must be a dihedral model" )
    end select

  end subroutine EmDee_add_dihedral

  !===================================================================================================

    subroutine EmDee_download( md, option, address ) bind(C,name="EmDee_download")
      type(tEmDee),      value      :: md
      character(c_char), intent(in) :: option(*)
      type(c_ptr),       value      :: address

      real(rb), pointer :: scalar, matrix(:,:)
      type(tData), pointer :: me
      character(sl) :: item

      call c_f_pointer( md%data, me )
      item = string(option)
      if (.not.c_associated(address)) call error( "download", "provided address is invalid" )

      select case (item)

        case ("box")
          call c_f_pointer( address, scalar )
          scalar = me%Lbox

        case ("coordinates")
          if (.not.associated( me%R )) call error( "download", "coordinates have not been allocated" )
          call c_f_pointer( address, matrix, [3,me%natoms] )
          !$omp parallel num_threads(me%nthreads)
          call download( omp_get_thread_num() + 1, me%R, matrix )
          !$omp end parallel

        case ("momenta")
          call c_f_pointer( address, matrix, [3,me%natoms] )
          !$omp parallel num_threads(me%nthreads)
          call get_momenta( omp_get_thread_num() + 1, matrix )
          !$omp end parallel

        case ("forces")
          call c_f_pointer( address, matrix, [3,me%natoms] )
          if (.not.me%forcesUpToDate(me%layer)) call EmDee_compute_forces( md )
!          if (me%kspace_active) then
!            call me % kspace % discount_rigid_pairs( me%layer, me%Lbox*[1,1,1], me%R, Forces = me%F )
!          end if
          !$omp parallel num_threads(me%nthreads)
          call download( omp_get_thread_num() + 1, me%F, matrix )
          !$omp end parallel

if (PRINT) print*, "downloaded from layer ", me%layer

        case ("centersOfMass")
          call c_f_pointer( address, matrix, [3,me%nbodies+me%nfree] )
          !$omp parallel num_threads(me%nthreads)
          call get_centers_of_mass( omp_get_thread_num() + 1 )
          !$omp end parallel

        case ("quaternions", "quatmom", "quattau")
          call c_f_pointer( address, matrix, [4,me%nbodies] )
          !$omp parallel num_threads(me%nthreads)
          call get_quaternions( omp_get_thread_num() + 1, item, matrix )
          !$omp end parallel

        case ("angmom", "bodycoord", "bodymom", "bodyforces", "torques", "inertia")
          call c_f_pointer( address, matrix, [3,me%nbodies] )
          !$omp parallel num_threads(me%nthreads)
          call get_body_properties( omp_get_thread_num() + 1, item, matrix )
          !$omp end parallel

        case default
          call error( "download", "invalid option" )

      end select

      contains
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        subroutine get_momenta( thread, Pext )
          integer,  intent(in)    :: thread
          real(rb), intent(inout) :: Pext(3,me%natoms)
          integer :: i
          forall (i = (thread - 1)*me%threadFreeAtoms + 1 : min(thread*me%threadFreeAtoms, me%nfree))
            Pext(:,me%free(i)) = me%P(:,me%free(i))
          end forall
          forall(i = (thread - 1)*me%threadBodies + 1 : min(thread*me%threadBodies,me%nbodies))
            Pext(:,me%body(i)%index) = me%body(i) % particle_momenta()
          end forall
        end subroutine get_momenta
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        subroutine download( thread, origin, destination )
          integer,  intent(in)    :: thread
          real(rb), intent(in)    :: origin(3,me%natoms)
          real(rb), intent(inout) :: destination(3,me%natoms)
          integer :: first, last
          first = (thread - 1)*me%threadAtoms + 1
          last = min(thread*me%threadAtoms, me%natoms)
          destination(:,first:last) = origin(:,first:last)
        end subroutine download
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        subroutine get_centers_of_mass( thread )
          integer,  intent(in) :: thread
          integer :: i
          forall (i = (thread - 1)*me%threadBodies + 1 : min(thread*me%threadBodies,me%nbodies) )
            matrix(:,i) = me%body(i)%rcm
          end forall
          forall(i = (thread - 1)*me%threadFreeAtoms + 1 : min(thread*me%threadFreeAtoms,me%nfree) )
            matrix(:,i+me%nbodies) = me%R(:,me%free(i))
          end forall
        end subroutine get_centers_of_mass
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        subroutine get_quaternions( thread, item, matrix )
          integer,       intent(in)  :: thread
          character(sl), intent(in)  :: item
          real(rb),      intent(out) :: matrix(:,:)

          integer :: i
          do i = (thread - 1)*me%threadBodies + 1, min(thread*me%threadBodies,me%nbodies)
            select case (item)
              case ("quaternions")
                matrix(:,i) = me%body(i)%q
              case ("quatmom")
                matrix(:,i) = me%body(i)%pi
              case ("quattau")
                matrix(:,i) = matmul(matrix_C(me%body(i)%q), two*me%body(i)%tau)
            end select
          end do
        end subroutine get_quaternions
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        subroutine get_body_properties( thread, item, matrix )
          integer,       intent(in)  :: thread
          character(sl), intent(in)  :: item
          real(rb),      intent(out) :: matrix(:,:)

          integer :: i
          do i = (thread - 1)*me%threadBodies + 1, min(thread*me%threadBodies,me%nbodies)
            select case (item)
              case ("angmom")
                matrix(:,i) = me%body(i)%omega
              case ("bodycoord")
                matrix(:,i) = me%body(i)%Rcm
              case ("bodymom")
                matrix(:,i) = me%body(i)%Pcm
              case ("bodyforces")
                matrix(:,i) = me%body(i)%F
              case ("torques")
                matrix(:,i) = me%body(i)%tau
              case ("inertia")
                matrix(:,i) = me%body(i)%MoI
              end select
          end do
        end subroutine get_body_properties
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    end subroutine EmDee_download

!===================================================================================================

  subroutine EmDee_upload( md, option, address ) bind(C,name="EmDee_upload")
    type(tEmDee),      intent(inout) :: md
    character(c_char), intent(in)    :: option(*)
    type(c_ptr),       value         :: address

    real(rb), allocatable :: twoKEt(:,:), twoKEr(:,:)
    real(rb),     pointer :: scalar, Vector(:), Matrix(:,:)
    type(tData),  pointer :: me
    character(sl) :: item

    call c_f_pointer( md%data, me )
    item = string(option)
    if (.not.c_associated(address)) call error( "upload", "provided address is invalid" )

    select case (item)

      case ("box")
        call c_f_pointer( address, scalar )
        if (.not.associated(me%Lbox)) allocate( me%Lbox )
        me%Lbox = scalar
        if (me%initialized) then
          me%forcesUpToDate = .false.
          me%layerEnergy%UpToDate = .false.
          md%Energy%UpToDate = .false.
        else if (associated( me%R )) then
          call initialize_system()
        end if

      case ("coordinates")
        if (.not.associated( me%R )) allocate( me%R(3,me%natoms) )
        call c_f_pointer( address, Matrix, [3,me%natoms] )
        !$omp parallel num_threads(me%nthreads)
        call upload( omp_get_thread_num() + 1, Matrix, me%R )
        !$omp end parallel
        if (me%initialized) then
          me%forcesUpToDate = .false.
          me%layerEnergy%UpToDate = .false.
          md%Energy%UpToDate = .false.
          if (md%Options%AutoBodyUpdate) then
            !$omp parallel num_threads(me%nthreads)
            call update_rigid_bodies( me, omp_get_thread_num() + 1 )
            !$omp end parallel
          end if
        else if (associated(me%Lbox)) then
          call initialize_system()
        end if

      case ("momenta")
        if (.not.me%initialized) call error( "upload", "box and coordinates have not been defined" )
        call c_f_pointer( address, Matrix, [3,me%natoms] )
        allocate( twoKEt(3,me%nthreads), twoKEr(3,me%nthreads) )
        !$omp parallel num_threads(me%nthreads)
        block
          integer :: thread
          thread = omp_get_thread_num() + 1
          call assign_momenta( me, thread, Matrix, twoKEt(:,thread), twoKEr(:,thread) )
        end block
        !$omp end parallel
        md%Kinetic%RotPart = half*sum(TwoKEr,2)
        md%Kinetic%TransPart = half*sum(twoKEt,2)
        md%Kinetic%Rotational = sum(md%Kinetic%RotPart)
        md%Kinetic%Total = sum(md%Kinetic%TransPart) + md%Kinetic%Rotational
        md%Kinetic%ShadowKinetic = md%Kinetic%Total
        md%Kinetic%ShadowRotational = md%Kinetic%Rotational
        md%Kinetic%UpToDate = .true.

      case ("forces")
        if (.not.me%initialized) call error( "upload", "box and coordinates have not been defined" )
        call c_f_pointer( address, Matrix, [3,me%natoms] )
        !$omp parallel num_threads(me%nthreads)
        call upload( omp_get_thread_num() + 1, Matrix, me%F )
        !$omp end parallel

if (PRINT) print*, "uploaded to layer ", me%layer

      case ("charges")
        if (me%initialized) &
          call error( "upload", "cannot set charges after box and coordinates initialization" )
        call c_f_pointer( address, Vector, [me%natoms] )
        !$omp parallel num_threads(me%nthreads)
        call assign_charges( omp_get_thread_num() + 1, Vector )
        !$omp end parallel
        me%forcesUpToDate = .false.
        me%layerEnergy%UpToDate = .false.
        md%Energy%UpToDate = .false.

      case default
        call error( "upload", "invalid option" )

    end select

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine assign_charges( thread, Qext )
        integer,  intent(in) :: thread
        real(rb), intent(in) :: Qext(me%natoms)
        integer :: first, last
        first = (thread - 1)*me%threadAtoms + 1
        last = min(thread*me%threadAtoms, me%natoms)
        me%charge(first:last) = Qext(first:last)
        me%charged(first:last) = abs(Qext(first:last)) > epsilon(one)
      end subroutine assign_charges
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine upload( thread, origin, destination )
        integer,  intent(in)    :: thread
        real(rb), intent(in)    :: origin(3,me%natoms)
        real(rb), intent(inout) :: destination(3,me%natoms)
        integer :: first, last
        first = (thread - 1)*me%threadAtoms + 1
        last = min(thread*me%threadAtoms, me%natoms)
        destination(:,first:last) = origin(:,first:last)
      end subroutine upload
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine initialize_system()
        integer :: layer
        call perform_initialization( me, md%DOF, md%RotDoF )
        do layer = me%nlayers, 1, -1
          call EmDee_switch_model_layer( md, layer )
          call EmDee_compute_forces( md )
        end do
      end subroutine initialize_system
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine EmDee_upload

  !===================================================================================================

    subroutine EmDee_switch_model_layer( md, layer ) bind(C,name="EmDee_switch_model_layer")
      type(tEmDee), intent(inout) :: md
      integer(ib),  value         :: layer

      type(tData), pointer :: me

      call c_f_pointer( md%data, me )
      if (layer /= me%layer) then
        if ((layer < 1).or.(layer > me%nlayers)) then
          call error( "model layer switch", "selected layer is out of range" )
        end if
        me%layer = layer
        md%Energy = me%layerEnergy(layer)
        md%Virial = me%layerVirial(layer)
        me%F => me%layerF(:,:,layer)
if (PRINT) print*, "switched to layer ", me%layer
      end if

    end subroutine EmDee_switch_model_layer

!===================================================================================================

  subroutine EmDee_random_momenta( md, kT, adjust, seed ) bind(C,name="EmDee_random_momenta")
    type(tEmDee), intent(inout) :: md
    real(rb),     value         :: kT
    logical(lb),  value         :: adjust
    integer(ib),  value         :: seed

    integer  :: i, j
    real(rb) :: twoKEt(3), TwoKEr(3)
    type(tData), pointer :: me

    call c_f_pointer( md%data, me )
    if (me%random%seeding_required) call me % random % setup( seed )
    twoKEt = zero
    TwoKEr = zero
    associate (rng => me%random)
      if (me%nbodies /= 0) then
        if (.not.me%initialized) call error( "random_momenta", "coordinates have not defined" )
        do i = 1, me%nbodies
          associate (b => me%body(i))
            b%pcm = sqrt(b%mass*kT)*[rng%normal(), rng%normal(), rng%normal()]
            call b%assign_momenta( sqrt(b%invMoI*kT)*[rng%normal(), rng%normal(), rng%normal()] )
            twoKEt = twoKEt + b%invMass*b%pcm**2
            TwoKEr = TwoKEr + b%MoI*b%omega**2
          end associate
        end do
      end if
      do j = 1, me%nfree
        i = me%free(j)
        me%P(:,i) = sqrt(me%mass(i)*kT)*[rng%normal(), rng%normal(), rng%normal()]
        twoKEt = twoKEt + me%invMass(i)*me%P(:,i)**2
      end do
    end associate
    if (adjust) call adjust_momenta

    md%Kinetic%RotPart = half*TwoKEr
    md%Kinetic%TransPart = half*twoKEt
    md%Kinetic%Rotational = sum(md%Kinetic%RotPart)
    md%Kinetic%Total = sum(md%Kinetic%TransPart) + md%Kinetic%Rotational
    md%Kinetic%ShadowKinetic = md%Kinetic%Total
    md%Kinetic%ShadowRotational = md%Kinetic%Rotational
    md%Kinetic%UpToDate = .true.

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine adjust_momenta
        integer  :: i, j
        real(rb) :: vcm(3), factor
        forall (i=1:3) vcm(i) = (sum(me%P(i,me%free)) + sum(me%body%pcm(i)))/me%totalMass
        twoKEt = zero
        do j = 1, me%nfree
          i = me%free(j)
          me%P(:,i) = me%P(:,i) - me%mass(i)*vcm
          twoKEt = twoKEt + me%invMass(i)*me%P(:,i)**2
        end do
        do i = 1, me%nbodies
          me%body(i)%pcm = me%body(i)%pcm - me%body(i)%mass*vcm
          twoKEt = twoKEt + me%body(i)%invMass*me%body(i)%pcm**2
        end do
        factor = sqrt((3*me%nfree + sum(me%body%dof) - 3)*kT/sum(twoKEt + TwoKEr))
        me%P(:,me%free) = factor*me%P(:,me%free)
        do i = 1, me%nbodies
          associate( b => me%body(i) )
            b%pcm = factor*b%pcm
            call b%assign_momenta( factor*b%omega )
          end associate
        end do
        twoKEt = factor*factor*twoKEt
        TwoKEr = factor*factor*TwoKEr
      end subroutine adjust_momenta
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine EmDee_random_momenta

!===================================================================================================

  subroutine EmDee_boost( md, lambda, alpha, dt ) bind(C,name="EmDee_boost")
    type(tEmDee), intent(inout) :: md
    real(rb),     value         :: lambda, alpha, dt

    real(rb) :: CP, CF
    real(rb), allocatable :: twoKEt(:,:), twoKEr(:,:)
    type(tData), pointer :: me

    call c_f_pointer( md%data, me )

    CF = phi(alpha*dt)*dt
    CP = one - alpha*CF
    CF = lambda*CF

    if ((lambda /= zero).and.(.not.me%forcesUpToDate(me%layer))) then
      call EmDee_compute_forces( md )
    end if

    if (md%Options%Compute) allocate( twoKEt(3,me%nthreads), twoKEr(3,me%nthreads) )
    !$omp parallel num_threads(me%nthreads)
    block
      integer :: thread
      thread = omp_get_thread_num() + 1
      call boost( me, thread, CP, CF, me%F, md%Options%translate, md%Options%rotate )
      if (md%Options%Compute) then
        call kinetic_energies( me, thread, md%Options%translate, md%Options%rotate, &
                               twoKEt(:,thread), twoKEr(:,thread) )
      end if
    end block
    !$omp end parallel

    if (md%Options%Compute) then
      if (md%Options%translate) md%Kinetic%TransPart = half*sum(twoKEt,2)
      if (md%Options%rotate) then
        md%Kinetic%RotPart = half*sum(twoKEr,2)
        md%Kinetic%Rotational = sum(md%Kinetic%RotPart)
      end if
      md%Kinetic%Total = sum(md%Kinetic%TransPart) + md%Kinetic%Rotational
    end if
    md%Kinetic%UpToDate = md%Options%Compute

if (PRINT) print*, "BOOST - layer = ", me%layer, " | dt = ", dt

  end subroutine EmDee_boost

!===================================================================================================

  subroutine EmDee_displace( md, lambda, alpha, dt ) bind(C,name="EmDee_displace")
    type(tEmDee), intent(inout) :: md
    real(rb),     value         :: lambda, alpha, dt

    real(rb) :: CR, CP
    type(tData), pointer :: me

    call c_f_pointer( md%data, me )

    md%Time%Motion = md%Time%Motion - omp_get_wtime()

    if (alpha /= zero) then
      CP = phi(alpha*dt)*dt
      CR = one - alpha*CP
      me%Lbox = cR*me%Lbox
    else
      CP = dt
      CR = one
    end if
    CP = lambda*CP

    associate ( Opt => md%Options )
      !$omp parallel num_threads(me%nthreads)
      call move( me, omp_get_thread_num() + 1, CR, CP, dt, &
                                  Opt%translate, Opt%rotate, Opt%rotationMode )
      !$omp end parallel
    end associate

    me%forcesUpToDate = .false.
    me%layerEnergy%UpToDate = .false.
    md%Energy%UpToDate = .false.

    md%Time%Motion = md%Time%Motion + omp_get_wtime()

if (PRINT) print*, "MOVE dt = ", dt

  end subroutine EmDee_displace

!===================================================================================================

  subroutine EmDee_verlet_step( md, dt ) bind(C,name="EmDee_verlet_step")
    type(tEmDee), intent(inout) :: md
    real(rb),     value         :: dt

    real(rb) :: Us, Ks_t, Ks_r, dt_2
    real(rb), allocatable :: twoKEt(:,:), twoKEr(:,:)
    real(rb), allocatable :: r0(:,:), q0(:,:), s0(:,:)
    type(tData),  pointer :: me

    dt_2 = half*dt
    call c_f_pointer( md%data, me )

    if (md%Options%Compute) then
      Us = zero
      Ks_t = zero
      Ks_r = zero
      allocate( twoKEt(3,me%nthreads), twoKEr(3,me%nthreads), &
                r0(3,me%nbodies), q0(4,me%nbodies), s0(3,me%nfree) )
    end if

    !$omp parallel num_threads(me%nthreads)
    call pre_force(omp_get_thread_num() + 1)
    !$omp end parallel

    call EmDee_compute_forces( md )

    !$omp parallel num_threads(me%nthreads) reduction(+:Us,Ks_t,Ks_r)
    call post_force(omp_get_thread_num() + 1, Us, Ks_t, Ks_r)
    !$omp end parallel

    if (md%Options%Compute) then
      md%Kinetic%TransPart = half*sum(twoKEt,2)
      md%Kinetic%RotPart = half*sum(twoKEr,2)
      md%Kinetic%Rotational = sum(md%Kinetic%RotPart)
      md%Kinetic%Total = sum(md%Kinetic%TransPart) + md%Kinetic%Rotational
      md%Energy%ShadowPotential = md%Energy%ShadowPotential - dt*dt*Us/24.0_rb
      md%Kinetic%ShadowRotational = Ks_r/(6.0_rb*dt)
      md%Kinetic%ShadowKinetic = (Ks_t + Ks_r)/(6.0_rb*dt)
    end if
    md%Kinetic%UpToDate = md%Options%Compute

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      function virtual_rotation( b, dt ) result ( q )
        type(tBody), intent(in)  :: b
        real(rb),    intent(in)  :: dt
        real(rb)                 :: q(4)
        type(tBody) :: c
        c = b
        call c % assign_momenta( c%pi + matmul(matrix_C(c%q), dt*c%tau) )
        if (md%Options%RotationMode /= 0) then
          call c % rotate_no_squish( dt, md%Options%RotationMode )
        else
          call c % rotate_exact( dt )
        end if
        q = c%q
      end function virtual_rotation
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine pre_force( thread )
        integer, intent(in) :: thread
        integer :: i, j
        if (md%Options%Compute) then
          do i = (thread - 1)*me%threadBodies + 1, min(thread*me%threadBodies,me%nbodies)
            associate (b => me%body(i))
              r0(:,i) = 2.5_rb*b%rcm + dt_2*b%invMass*(b%pcm - dt_2*b%F)
              q0(:,i) = half*virtual_rotation( b, -dt ) - three*b%q
            end associate
          end do
          do i = (thread - 1)*me%threadFreeAtoms + 1, min(thread*me%threadFreeAtoms, me%nfree)
            j = me%free(i)
            s0(:,i) = 2.5_rb*me%R(:,j) + dt_2*me%invMass(j)*(me%P(:,j) - dt_2*me%F(:,j))
          end do
        end if
        call boost( me, thread, one, dt_2, me%F, TRUE, TRUE )
        call move( me, thread, one, dt, dt, TRUE, TRUE, md%Options%rotationMode )
      end subroutine pre_force
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine post_force( thread, Us, Ks_t, Ks_r )
        integer,  intent(in)    :: thread
        real(rb), intent(inout) :: Us, Ks_t, Ks_r
        integer :: i, j
        real(rb) :: tau_b(3), rdot(3), qdot(4)
        call boost( me, thread, one, dt_2, me%F, TRUE, TRUE )
        if (md%Options%Compute) then
          call kinetic_energies( me, thread, TRUE, TRUE, twoKEt(:,thread), twoKEr(:,thread) )
          do i = (thread - 1)*me%threadBodies + 1, min(thread*me%threadBodies,me%nbodies)
            associate (b => me%body(i))
              rdot =  2.5_rb*b%rcm + dt*b%invMass*(b%pcm + dt_2*b%F) - r0(:,i)
              Ks_t = Ks_t + sum(rdot*b%pcm)
              qdot = q0(:,i) + 1.5_rb*b%q + virtual_rotation(b, dt)
              Ks_r = Ks_r + sum((qdot - sum(qdot*b%q)*b%q)*b%pi)
              tau_b = matmul(matmul(matrix_Bt(b%q), matrix_C(b%q)), b%tau)
              Us = Us + b%invMass*sum(b%F**2) + sum(b%invMoI*tau_b**2)
            end associate
          end do
          do i = (thread - 1)*me%threadFreeAtoms + 1, min(thread*me%threadFreeAtoms, me%nfree)
            j = me%free(i)
            rdot =  2.5_rb*me%R(:,j) + dt*me%invMass(j)*(me%P(:,j) + dt_2*me%F(:,j)) - s0(:,i)
            Ks_t = Ks_t + sum(rdot*me%P(:,j))
            Us = Us + me%invMass(j)*sum(me%F(:,j)**2)
          end do
        end if
      end subroutine post_force
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine EmDee_verlet_step

!===================================================================================================

  subroutine EmDee_compute_forces( md ) bind(C,name="EmDee_compute_forces")
    type(tEmDee), intent(inout) :: md

    integer, parameter :: pair=1, coul=2, long=3, bond=4, angle=5
    real(rb) :: E(5), W(5)
    real(rb) :: time
    logical(lb) :: compute
    real(rb), allocatable :: Rs(:,:), Fs(:,:,:)
    type(tData),  pointer :: me

    call c_f_pointer( md%data, me )

    allocate( Rs(3,me%natoms), Fs(3,me%natoms,me%nthreads) )
    Rs = me%R/me%Lbox

    call handle_neighbor_lists( me, md%builds, md%Time%Neighbor, Rs )

    md%Time%Pair = md%Time%Pair - omp_get_wtime()
    compute = md%Options%Compute
    E = zero
    W = zero
    !$omp parallel num_threads(me%nthreads) reduction(+:E,W)
    block
      integer :: thread
      thread = omp_get_thread_num() + 1
      associate( F => Fs(:,:,thread) )
        call compute_pairs( me, thread, compute, Rs, F, E(pair), E(coul), W(pair), W(coul) )
        call compute_bonds( me, thread, Rs, F, E(bond), W(bond), E(coul) )
        call compute_angles( me, thread, Rs, F, E(angle), W(angle), E(coul) )
      end associate
    end block
    !$omp end parallel
    me%F = sum(Fs,3)

    if (me%coul(me%layer)%model%requires_kspace) then
      call compute_kspace( me, Rs, E(long), me%F )
      W(long) = E(coul) + E(long) - W(coul)
    end if

    md%Virial%Total = sum(W)
    if (me%nbodies /= 0) then
      md%Virial%Body = rigid_body_virial( me )
      md%Virial%Total = md%Virial%Total + md%Virial%Body
    end if

    if (compute) then
      md%Energy%Dispersion = E(pair)
      md%Energy%Coulomb = E(coul) + E(long)
      md%Energy%Bond = E(bond)
      md%Energy%Angle = E(angle)
      md%Energy%Potential = sum(E)
      md%Energy%ShadowPotential = md%Energy%Potential
    end if

    md%Energy%UpToDate = compute
    me%forcesUpToDate(me%layer) = .true.
    me%layerEnergy(me%layer) = md%Energy
    me%layerVirial(me%layer) = md%Virial
    time = omp_get_wtime()
    md%Time%Pair = md%Time%Pair + time
    md%Time%Total = time - me%startTime

if (PRINT) print*, "FORCE - layer = ", me%layer

  end subroutine EmDee_compute_forces

!===================================================================================================

  subroutine EmDee_rdf( md, bins, pairs, itype, jtype, g ) bind(C,name="EmDee_rdf")
    type(tEmDee), value       :: md
    integer(ib),  value       :: pairs, bins
    integer(ib),  intent(in)  :: itype(pairs), jtype(pairs)
    real(rb),     intent(out) :: g(bins,pairs)

    character(*), parameter :: task = "radial distribution calculation"
    real(rb),     parameter :: Pi4_3 = 4.188790204786391_rb

    integer :: maxtype, i, j, bin, pair
    real(rb) :: invL, invL2

    logical,  allocatable :: hasPair(:), pairOn(:,:)
    integer,  allocatable :: pairCount(:,:,:), N(:)
    real(rb), allocatable :: Rs(:,:), rdf(:,:)
    type(tData),  pointer :: me

    call c_f_pointer( md%data, me )
    if (.not.ranged( [itype,jtype], me%ntypes )) then
      call error( task, "at least one provided type index is out of range" )
    end if

    allocate( hasPair(me%ntypes), source = .false. )
    hasPair(itype) = .true.
    hasPair(jtype) = .true.

    allocate( pairOn(me%ntypes,me%ntypes), source = .false. )
    pairOn(itype,jtype) = .true.
    pairOn(jtype,itype) = .true.

    maxtype = maxval([itype,jtype])
    allocate( pairCount(bins,symm1D(maxtype,maxtype),me%nthreads), Rs(3,me%natoms) )

    invL = one/me%Lbox
    invL2 = invL*invL

    !$omp parallel num_threads(me%nthreads)
    block
      integer :: thread, a1, aN
      thread = omp_get_thread_num() + 1
      a1 = (thread - 1)*me%threadAtoms + 1
      aN = min(thread*me%threadAtoms, me%natoms)
      Rs(:,a1:aN) = invL*me%R(:,a1:aN)
      !$omp barrier
      call count_pairs( thread, Rs )
    end block
    !$omp end parallel

    rdf = sum(pairCount,3)/(Pi4_3*(me%Rc*invL/bins)**3)
    forall (bin=1:bins) rdf(bin,:) = rdf(bin,:)/(3*bin*(bin - 1) + 1)
    allocate( N(maxtype) )
    forall(i=1:maxtype, hasPair(i)) N(i) = count(me%atomType == i)
    do pair = 1, pairs
      i = itype(pair)
      j = jtype(pair)
      if (i == j) then
        g(:,pair) = two*rdf(:,symm1D(i,j))/(N(i)*N(j))
      else
        g(:,pair) = rdf(:,symm1D(i,j))/(N(i)*N(j))
      end if
    end do

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine count_pairs( thread, Rs )
        integer,     intent(in) :: thread
        real(rb),    intent(in) :: Rs(3,me%natoms)

        integer  :: firstAtom, lastAtom, k, i, itype, m, j, jtype, pair, bin
        real(rb) :: Rc2, binsByRc, Ri(3), Rij(3), r2

        Rc2 = me%RcSq*invL2
        binsByRc = bins/(me%Rc*invL)
        pairCount(:,:,thread) = 0
        associate ( neighbor => me%neighbor(thread) )
          firstAtom = me%cellAtom%first(me%threadCell%first(thread))
          lastAtom = me%cellAtom%last(me%threadCell%last(thread))
          do k = firstAtom, lastAtom
            i = me%cellAtom%item(k)
            itype = me%atomType(i)
            if (hasPair(itype)) then
              Ri = Rs(:,i)
              do m = neighbor%first(i), neighbor%last(i)
                j = neighbor%item(m)
                jtype = me%atomType(j)
                if (pairOn(itype,jtype)) then
                  Rij = pbc(Ri - Rs(:,j))
                  r2 = sum(Rij*Rij)
                  if (r2 < Rc2) then
                    pair = symm1D(itype,jtype)
                    bin = int(sqrt(r2)*binsByRc) + 1
                    pairCount(bin,pair,thread) = pairCount(bin,pair,thread) + 1
                  end if
                end if
              end do
            end if
          end do
        end associate

      end subroutine count_pairs
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      elemental function pbc( x )
        real(rb), intent(in) :: x
        real(rb)              :: pbc
        pbc = x - anint(x)
      end function pbc
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine EmDee_rdf

!===================================================================================================
!                              A U X I L I A R Y   P R O C E D U R E S
!===================================================================================================
  !
  ! subroutine update_forces( md, layer )
  !   type(tEmDee), intent(inout) :: md
  !   integer,      intent(in)    :: layer
  !
  !   real(rb) :: DEpair, DEcoul, DWpair, DWcoul, time
  !   real(rb), allocatable :: Rs(:,:), DF(:,:,:)
  !   type(tData),  pointer :: me
  !
  !   call c_f_pointer( md%data, me )
  !   md%Time%Pair = md%Time%Pair - omp_get_wtime()
  !
  !   allocate( Rs(3,me%natoms), DF(3,me%natoms,me%nthreads) )
  !   Rs = (one/me%Lbox)*me%R
  !
  !   !$omp parallel num_threads(me%nthreads) reduction(+:DEpair,DEcoul,DWpair,DWcoul)
  !   block
  !     integer :: thread
  !     thread = omp_get_thread_num() + 1
  !     call update_pairs( me, thread, Rs, DF(:,:,thread), DEpair, DEcoul, DWpair, DWcoul, layer )
  !   end block
  !   !$omp end parallel
  !
  !   me%F = me%F + sum(DF,3)
  !   md%Energy%Potential = md%Energy%Potential + DEpair
  !   md%Energy%Coulomb = md%Energy%Coulomb + DEcoul
  !   md%Virial = md%Virial + DWpair + DWcoul
  !
  !   time = omp_get_wtime()
  !   md%Time%Pair = md%Time%Pair + time
  !   md%Time%Total = time - me%startTime
  !
  ! end subroutine update_forces
  !
!===================================================================================================

end module EmDeeCode
