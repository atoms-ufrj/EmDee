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

! TODO: Include kspace terms in subroutine update_forces.
! TODO: Discount long-range electrostatic terms of bonds, angles, and dihedrals
! TODO: Change name of exported Julia wrapper
! TODO: Replace type(c_ptr) arguments by array arguments when possible
! TODO: Create indexing for having sequential body particles and free particles in arrays

module EmDeeCode

use EmDeeData
use neighbor_lists

implicit none

private

character(11), parameter :: VERSION = "03 Mar 2017"

type, bind(C), public :: tOpts
  logical(lb) :: Translate            ! Flag to activate/deactivate translations
  logical(lb) :: Rotate               ! Flag to activate/deactivate rotations
  integer(ib) :: RotationMode         ! Algorithm used for free rotation of rigid bodies
  real(rb)    :: Lambda_R             ! Momentum-multiplying constant in position equations
  real(rb)    :: Alpha_R              ! Position-multiplying constant in position equations
  real(rb)    :: Lambda_P             ! Force-multiplying constant in momentum equations
  real(rb)    :: Alpha_P              ! Momentum-multiplying constant in momentum equations
end type tOpts

type, bind(C), public :: tEnergy
  real(rb)    :: Potential            ! Total potential energy of the system
  real(rb)    :: Dispersion           ! Dispersion (vdW) part of the potential energy
  real(rb)    :: Coulomb              ! Electrostatic part of the potential energy
  real(rb)    :: Fourier              ! Reciprocal part of the electrostatic potential
  real(rb)    :: Kinetic              ! Total kinetic energy of the system
  real(rb)    :: TransPart(3)         ! Translational kinetic energy at each dimension
  real(rb)    :: Rotational           ! Total rotational kinetic energy of the system
  real(rb)    :: RotPart(3)           ! Rotational kinetic energy around each principal axis
  type(c_ptr) :: LayerPotential       ! Vector with multilayer potential energy components
  type(c_ptr) :: LayerDispersion      ! Vector with multilayer dispersion energy components
  type(c_ptr) :: LayerCoulomb         ! Vector with multilayer coulombic energy components
  type(c_ptr) :: LayerFourier         ! Vector with multilayer reciprocal energy components
  logical(lb) :: Compute              ! Flag to activate/deactivate energy computations
  logical(lb) :: UpToDate             ! Flag to attest whether energies have been computed
end type tEnergy

type, bind(C), public :: tTime
  real(rb) :: Pair
  real(rb) :: FastPair
  real(rb) :: Neighbor
  real(rb) :: Total
end type tTime

type, bind(C), public :: tEmDee
  integer(ib)   :: Builds             ! Number of neighbor list builds
  type(tTime)   :: Time
  type(tEnergy) :: Energy             ! All energy terms
  real(rb)      :: Virial             ! Total internal virial of the system
  real(rb)      :: BodyVirial         ! Rigid body contribution to the internal virial
  integer(ib)   :: DoF                ! Total number of degrees of freedom
  integer(ib)   :: RotDoF             ! Number of rotational degrees of freedom
  type(c_ptr)   :: Data               ! Pointer to system data
  type(tOpts)   :: Options            ! List of options to change EmDee's behavior
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

    integer :: i
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
    me%skinSq = skin*skin
    me%natoms = N
    me%threadAtoms = (N + threads - 1)/threads
    me%kspace_active = .false.
    me%respa_active = .false.
    me%xRespaRcSq = zero

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
    allocate( me%P(3,N), me%F(3,N), me%R0(3,N), source = zero )
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

    ! Allocate variables related to model layers:
    me%layer = 1
    me%other_layer = [(i,i=2,me%nlayers)]
    allocate( me%threadEpair(layers,threads), source = zero )
    allocate( me%threadEcoul(layers,threads), source = zero )
    allocate( me%threadElong(layers,threads), source = zero )
    allocate( me%Epair(layers), source = zero )
    allocate( me%Ecoul(layers), source = zero )
    allocate( me%Elong(layers), source = zero )
    allocate( me%Energy(layers), source = zero )

    ! Set up mutable entities:
    EmDee_system % builds = 0

    EmDee_system % Time % Pair = zero
    EmDee_system % Time % FastPair = zero
    EmDee_system % Time % Neighbor = zero
    EmDee_system % Time % Total = zero

    EmDee_system % Energy % Potential = zero
    EmDee_system % Energy % Dispersion = zero
    EmDee_system % Energy % Coulomb = zero
    EmDee_system % Energy % Fourier = zero
    EmDee_system % Energy % Kinetic = zero
    EmDee_system % Energy % TransPart = zero
    EmDee_system % Energy % Rotational = zero
    EmDee_system % Energy % RotPart = zero
    EmDee_system % Energy % LayerPotential = c_loc(me%Energy(1))
    EmDee_system % Energy % LayerDispersion = c_loc(me%Epair(1))
    EmDee_system % Energy % LayerCoulomb = c_loc(me%Ecoul(1))
    EmDee_system % Energy % LayerFourier = c_loc(me%Elong(1))
    EmDee_system % Energy % Compute = .true.
    EmDee_system % Energy % UpToDate = .false.

    EmDee_system % Virial = zero
    EmDee_system % BodyVirial = zero
    EmDee_system % DoF = 3*(N - 1)
    EmDee_system % RotDoF = 0
    EmDee_system % data = c_loc(me)
    EmDee_system % Options % translate = .true.

    EmDee_system % Options % Rotate = .true.
    EmDee_system % Options % RotationMode = 0
    EmDee_system % Options % Lambda_R = one
    EmDee_system % Options % Alpha_R  = zero
    EmDee_system % Options % Lambda_P = one
    EmDee_system % Options % Alpha_P  = zero

  end function EmDee_system

!===================================================================================================

  subroutine EmDee_switch_model_layer( md, layer ) bind(C,name="EmDee_switch_model_layer")
    type(tEmDee), value :: md
    integer(ib),  value :: layer

    integer :: i
    type(tData), pointer :: me

    call c_f_pointer( md%data, me )
    if (layer /= me%layer) then
      if ((layer < 1).or.(layer > me%nlayers)) then
        call error( "model layer switch", "selected layer is out of range" )
      end if
      if (me%initialized) call update_forces( md, layer )
      me%layer = layer
      me%other_layer = [(i,i=1,layer-1), (i,i=layer+1,me%nlayers)]
    end if

  end subroutine EmDee_switch_model_layer

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

    write(C,'(I5)') me%nlayers
    C = adjustl(C)
    do layer = 1, me%nlayers
      if (.not.c_associated(model(layer))) then
        call error( task, trim(C)//" valid pair models must be provided" )
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

    integer :: layer
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
          me%coul(layer) = container
          call me % coul(layer) % model % shifting_setup( me%Rc )
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

    integer :: layer
    character(5) :: C
    type(tData), pointer :: me
    type(modelContainer), pointer :: container

    call c_f_pointer( md%data, me )
    if (me%initialized) call error( task, "cannot set model after coordinates have been defined" )

    write(C,'(I5)') me%nlayers
    C = adjustl(C)
    do layer = 1, me%nlayers
      if (.not.c_associated(model(layer))) then
        call error( task, trim(C)//" valid coulomb models must be provided" )
      end if
      call c_f_pointer( model(layer), container )
      select type ( coul => container%model )
        class is (cCoulModel)
          me%coul(layer) = container
          call me % coul(layer) % model % shifting_setup( me%Rc )
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

  subroutine EmDee_set_respa( md, Rc, Npair, Nbond ) bind(C,name="EmDee_set_respa")
    type(tEmDee), intent(inout) :: md
    real(rb),     value         :: Rc
    integer(ib),  value         :: Npair, Nbond

    type(tData), pointer :: me

    call c_f_pointer( md%data, me )

    if (me%initialized) call error( "set RESPA", "the system has already been initialized" )
    if (Rc >= me%Rc) call error( "set RESPA", "provided internal cutoff is too large" )

    me%respa_active = .true.
    me%respaRc = Rc
    me%respaRcSq = Rc**2
    me%xRespaRcSq = (Rc + me%skin)**2
    me%fshift = one/me%respaRcSq
    me%Npair = Npair
    me%Nbond = Nbond

  end subroutine EmDee_set_respa

!===================================================================================================

  subroutine EmDee_upload( md, option, address ) bind(C,name="EmDee_upload")
    type(tEmDee),      intent(inout) :: md
    character(c_char), intent(in)    :: option(*)
    type(c_ptr),       value         :: address

    real(rb), allocatable :: twoKEt(:,:), twoKEr(:,:)
    real(rb), pointer :: scalar, Vector(:), Matrix(:,:)
    type(tData), pointer :: me
    character(sl) :: item

    call c_f_pointer( md%data, me )
    item = string(option)
    if (.not.c_associated(address)) call error( "upload", "provided address is invalid" )

    select case (item)

      case ("box")
        call c_f_pointer( address, scalar )
        me%Lbox = scalar
        me%Lbox3 = scalar
        me%invL = one/scalar
        me%invL2 = me%invL**2
        if (.not.me%initialized) then
          me%initialized = allocated( me%R )
          if (me%initialized) call perform_initialization( me, md%DOF, md%RotDoF )
        end if
        if (me%initialized) call compute_all_forces( md )

      case ("coordinates")
        if (.not.allocated( me%R )) allocate( me%R(3,me%natoms) )
        call c_f_pointer( address, Matrix, [3,me%natoms] )
        !$omp parallel num_threads(me%nthreads)
        call upload( omp_get_thread_num() + 1, Matrix, me%R )
        !$omp end parallel
        if (.not.me%initialized) then
          me%initialized = me%Lbox > zero
          if (me%initialized) call perform_initialization( me, md%DOF, md%RotDoF )
        end if
        if (me%initialized) call compute_all_forces( md )

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
        md%Energy%RotPart = half*sum(TwoKEr,2)
        md%Energy%TransPart = half*sum(twoKEt,2)
        md%Energy%Rotational = sum(md%Energy%RotPart)
        md%Energy%Kinetic = sum(md%Energy%TransPart) + md%Energy%Rotational

      case ("forces")
        if (.not.me%initialized) call error( "upload", "box and coordinates have not been defined" )
        call c_f_pointer( address, Matrix, [3,me%natoms] )
        !$omp parallel num_threads(me%nthreads)
        call upload( omp_get_thread_num() + 1, Matrix, me%F )
        !$omp end parallel

      case ("charges")
        if (me%initialized) &
          call error( "upload", "cannot set charges after box and coordinates initialization" ) 
        call c_f_pointer( address, Vector, [me%natoms] )
        !$omp parallel num_threads(me%nthreads)
        call assign_charges( omp_get_thread_num() + 1, Vector )
        !$omp end parallel
        if (me%initialized) call compute_all_forces( md )

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
  end subroutine EmDee_upload

!===================================================================================================

  subroutine EmDee_download( md, option, address ) bind(C,name="EmDee_download")
    type(tEmDee),      value      :: md
    character(c_char), intent(in) :: option(*)
    type(c_ptr),       value      :: address

    real(rb), pointer :: scalar, vector(:), matrix(:,:)
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
        if (.not.allocated( me%R )) call error( "download", "coordinates have not been allocated" )
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
        if (me%kspace_active) then
          call me % kspace % discount_rigid_pairs( me%layer, me%Lbox3, me%R, Forces = me%F )
        end if
        !$omp parallel num_threads(me%nthreads)
        call download( omp_get_thread_num() + 1, me%F, matrix )
        !$omp end parallel

      case ("multienergy")
        if (.not.md%Energy%UpToDate) call error( "download", "layer energies are outdated" )
        call c_f_pointer( address, vector, [me%nlayers] )
        vector = me%Epair + me%Ecoul

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
  end subroutine EmDee_download

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

    md%Energy%RotPart = half*TwoKEr
    md%Energy%TransPart = half*twoKEt
    md%Energy%Rotational = sum(md%Energy%RotPart)
    md%Energy%Kinetic = sum(md%Energy%TransPart) + md%Energy%Rotational

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

    if (md%Energy%compute) allocate( twoKEt(3,me%nthreads), twoKEr(3,me%nthreads) )
    !$omp parallel num_threads(me%nthreads)
    block
      integer :: thread
      thread = omp_get_thread_num() + 1
      call boost( me, thread, CP, CF, me%F, md%Options%translate, md%Options%rotate )
      if (md%Energy%compute) then
        call kinetic_energies( me, thread, md%Options%translate, md%Options%rotate, &
                               twoKEt(:,thread), twoKEr(:,thread) )
      end if
    end block
    !$omp end parallel

    if (md%Energy%compute) then
      if (md%Options%translate) md%Energy%TransPart = half*sum(twoKEt,2)
      if (md%Options%rotate) then
        md%Energy%RotPart = half*sum(twoKEr,2)
        md%Energy%Rotational = sum(md%Energy%RotPart)
      end if
      md%Energy%Kinetic = sum(md%Energy%TransPart) + md%Energy%Rotational
    end if

  end subroutine EmDee_boost

!===================================================================================================

  subroutine EmDee_displace( md, lambda, alpha, dt ) bind(C,name="EmDee_displace")
    type(tEmDee), intent(inout) :: md
    real(rb),     value         :: lambda, alpha, dt

    real(rb) :: CR, CP
    type(tData), pointer :: me

    call c_f_pointer( md%data, me )

    if (alpha /= zero) then
      CP = phi(alpha*dt)*dt
      CR = one - alpha*CP
      me%Lbox = cR*me%Lbox
      me%Lbox3 = me%Lbox
      me%InvL = one/me%Lbox
      me%invL2 = me%invL*me%invL
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

    call compute_all_forces( md )

  end subroutine EmDee_displace

!===================================================================================================

  subroutine EmDee_advance( md, alpha_R, alpha_P, dt ) bind(C,name="EmDee_advance")
    type(tEmDee), intent(inout) :: md
    real(rb),     value         :: alpha_R, alpha_P, dt

    integer(ib) :: mode, step
    logical(lb) :: trans, rot, changeBox
    real(rb)    :: dt_2, smalldt, smalldt_2, CPex, CFex, CPin, CFin, CR, CP
    real(rb), allocatable :: F_slow(:,:)
    real(rb), allocatable :: twoKEt(:,:), twoKEr(:,:)
    type(tData),  pointer :: me

    call c_f_pointer( md%data, me )

    if (md%Energy%compute) allocate( twoKEt(3,me%nthreads), twoKEr(3,me%nthreads) )

    dt_2 = half*dt
    trans = md%Options%translate
    rot = md%Options%rotate
    mode = md%Options%rotationMode

    CFex = phi(alpha_P*dt_2)*dt_2
    CPex = one - alpha_P*CFex
    changeBox = alpha_R /= zero

    if (me%respa_active) then

      smalldt = dt/me%Npair
      smalldt_2 = half*smalldt

      CFin = phi(alpha_P*smalldt_2)*smalldt_2
      CPin = one - alpha_P*CFin

      CP = phi(alpha_R*smalldt)*smalldt
      CR = one - alpha_R*CP

      allocate( F_slow(3,me%natoms) )

      if (me%fast_required) call compute_fast_forces( md )

      !$omp parallel num_threads(me%nthreads)
      block
        integer :: thread, i
        thread = omp_get_thread_num() + 1
        forall (i = (thread - 1)*me%threadAtoms + 1 : min(thread*me%threadAtoms, me%natoms))
          F_slow(:,i) = me%F(:,i) - me%F_fast(:,i)
        end forall
        !$omp barrier
        call boost( me, thread, CPex, CFex, F_slow, trans, rot )
      end block
      !$omp end parallel

!print*, "slow = ", sum(abs(F_slow)), "fast = ", sum(abs(me%F_fast))

      do step = 1, me%Npair
        !$omp parallel num_threads(me%nthreads)
        block
          integer :: thread
          thread = omp_get_thread_num() + 1
          call boost( me, thread, CPin, CFin, me%F_fast, trans, rot )
          !$omp barrier
          call move( me, thread, CR, CP, smalldt, trans, rot, mode )
        end block
        !$omp end parallel

        if (changeBox) then
          me%Lbox = CR*me%Lbox
          me%InvL = one/me%Lbox
          me%invL2 = me%invL*me%invL
          me%Lbox3 = me%Lbox
        end if

        call compute_fast_forces( md )

        !$omp parallel num_threads(me%nthreads)
        call boost( me, omp_get_thread_num() + 1, CPin, CFin, me%F_fast, trans, rot )
        !$omp end parallel
      end do

      call compute_all_forces( md )

      !$omp parallel num_threads(me%nthreads)
      block
        integer :: thread, i
        thread = omp_get_thread_num() + 1
        forall (i = (thread - 1)*me%threadAtoms + 1 : min(thread*me%threadAtoms, me%natoms))
          F_slow(:,i) = me%F(:,i) - me%F_fast(:,i)
        end forall
        !$omp barrier
        call boost( me, thread, CPex, CFex, F_slow, trans, rot )
        if (md%Energy%compute) then
          call kinetic_energies( me, thread, trans, rot, twoKEt(:,thread), twoKEr(:,thread) )
        end if
      end block
      !$omp end parallel

    else

      CP = phi(alpha_R*dt)*dt
      CR = one - alpha_R*CP

      !$omp parallel num_threads(me%nthreads)
      block
        integer :: thread
        thread = omp_get_thread_num() + 1
        call boost( me, thread, CPex, CFex, me%F, trans, rot )
        !$omp barrier
        call move( me, thread, CR, CP, dt, trans, rot, mode )
      end block
      !$omp end parallel

      if (changeBox) then
        me%Lbox = CR*me%Lbox
        me%InvL = one/me%Lbox
        me%invL2 = me%invL*me%invL
        me%Lbox3 = me%Lbox
      end if

      call compute_all_forces( md )

      !$omp parallel num_threads(me%nthreads)
      block
        integer :: thread
        thread = omp_get_thread_num() + 1
        call boost( me, thread, CPex, CFex, me%F, trans, rot )
        if (md%Energy%compute) then
          call kinetic_energies( me, thread, trans, rot, twoKEt(:,thread), twoKEr(:,thread) )
        end if
      end block
      !$omp end parallel

    end if

    if (md%Energy%compute) then
      if (md%Options%translate) md%Energy%TransPart = half*sum(twoKEt,2)
      if (md%Options%rotate) then
        md%Energy%RotPart = half*sum(twoKEr,2)
        md%Energy%Rotational = sum(md%Energy%RotPart)
      end if
      md%Energy%Kinetic = sum(md%Energy%TransPart) + md%Energy%Rotational
    end if

  end subroutine EmDee_advance

!===================================================================================================
!                              A U X I L I A R Y   P R O C E D U R E S
!===================================================================================================

  subroutine compute_all_forces( md )
    type(tEmDee), intent(inout) :: md

    real(rb) :: time, Epair, Ecoul, Elong, Wpair, Wcoul
    logical(lb) :: compute
    real(rb), allocatable :: Rs(:,:), Fs(:,:,:)
    type(tData),  pointer :: me

    call c_f_pointer( md%data, me )

    allocate( Rs(3,me%natoms), Fs(3,me%natoms,me%nthreads) )
    Rs = me%invL*me%R

    call handle_neighbor_lists( me, md%builds, md%Time%Neighbor, Rs )

    md%Time%Pair = md%Time%Pair - omp_get_wtime()
    compute = md%Energy%Compute
    !$omp parallel num_threads(me%nthreads) reduction(+:Epair,Ecoul,Elong,Wpair,Wcoul)
    block
      integer :: thread
      thread = omp_get_thread_num() + 1
      associate( F => Fs(:,:,thread) )
        call compute_pairs( me, thread, compute, Rs, F, Epair, Ecoul, Wpair, Wcoul )
        Elong = zero
        me%threadElong(:,thread) = zero
        if (me%bonds%exist) call compute_bonds( me, thread, Rs, F, Epair, Wpair )
        if (me%angles%exist) call compute_angles( me, thread, Rs, F, Epair, Wpair )
        if (me%dihedrals%exist) call compute_dihedrals( me, thread, Rs, F, Epair, Wpair )
      end associate
    end block
    !$omp end parallel
    me%F = sum(Fs,3)

    if (me%kspace_active) then
      call compute_kspace( me, compute, Rs, Elong, me%F )
      md%Virial = Wpair - 3.0_rb*(Ecoul + Elong)
    else
      md%Virial = Wpair + Wcoul
    end if

    if (me%nbodies /= 0) then
      md%BodyVirial = rigid_body_virial( me )
      md%Virial = md%Virial + md%BodyVirial
    end if

    if (compute) then
      md%Energy%Dispersion = Epair
      md%Energy%Coulomb = Ecoul + Elong
      md%Energy%Fourier = Elong
      md%Energy%Potential = Epair + md%Energy%Coulomb
      me%Epair = sum(me%threadEpair,2)
      me%Elong = sum(me%threadElong,2)
      me%Ecoul = sum(me%threadEcoul,2) + me%Elong
      me%Energy = me%Epair + me%Ecoul
    end if

    me%fast_required = me%respa_active

    md%Energy%UpToDate = compute
    time = omp_get_wtime()
    md%Time%Pair = md%Time%Pair + time
    md%Time%Total = time - me%startTime

  end subroutine compute_all_forces

!===================================================================================================

  subroutine compute_fast_forces( md )
    type(tEmDee), intent(inout) :: md

    real(rb) :: time
    real(rb), allocatable :: Rs(:,:), Fs(:,:,:)
    type(tData),  pointer :: me

    call c_f_pointer( md%data, me )

    allocate( Rs(3,me%natoms), Fs(3,me%natoms,me%nthreads) )
    Rs = me%invL*me%R

    call handle_neighbor_lists( me, md%builds, md%Time%Neighbor, Rs )

    md%Time%FastPair = md%Time%FastPair - omp_get_wtime()

    !$omp parallel num_threads(me%nthreads)
    block
      integer :: thread
      thread = omp_get_thread_num() + 1
      call compute_short_range_forces( me, thread, Rs, Fs(:,:,thread) )
    end block
    !$omp end parallel
    me%F_fast = sum(Fs,3)

    me%fast_required = .false.

    time = omp_get_wtime()
    md%Time%FastPair = md%Time%FastPair + time
    md%Time%Total = time - me%startTime

  end subroutine compute_fast_forces

!===================================================================================================

  subroutine update_forces( md, layer )
    type(tEmDee), intent(inout) :: md
    integer,      intent(in)    :: layer

    real(rb) :: DEpair, DEcoul, DWpair, DWcoul, time
    real(rb), allocatable :: Rs(:,:), DF(:,:,:)
    type(tData),  pointer :: me

    call c_f_pointer( md%data, me )
    md%Time%Pair = md%Time%Pair - omp_get_wtime()

    allocate( Rs(3,me%natoms), DF(3,me%natoms,me%nthreads) )
    Rs = me%invL*me%R

    !$omp parallel num_threads(me%nthreads) reduction(+:DEpair,DEcoul,DWpair,DWcoul)
    block
      integer :: thread
      thread = omp_get_thread_num() + 1
      call update_pairs( me, thread, Rs, DF(:,:,thread), DEpair, DEcoul, DWpair, DWcoul, layer )
    end block
    !$omp end parallel

    me%F = me%F + sum(DF,3)
    md%Energy%Potential = md%Energy%Potential + DEpair
    md%Energy%Coulomb = md%Energy%Coulomb + DEcoul
    md%Virial = md%Virial + DWpair + DWcoul

    time = omp_get_wtime()
    md%Time%Pair = md%Time%Pair + time
    md%Time%Total = time - me%startTime

  end subroutine update_forces

!===================================================================================================

end module EmDeeCode
