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

! TODO: Include coulomb multimodels in multilayer energy and virial computations
! TODO: Change name of exported Julia wrapper
! TODO: Add option of automatic versus manual force calculations
! TODO: Replace type(c_ptr) arguments by array arguments when possible
! TODO: Create indexing for having sequential body particles and free particles in arrays

module EmDeeCode

use EmDeeData

implicit none

private

character(11), parameter :: VERSION = "03 Feb 2017"

type, bind(C) :: tOpts
  logical(lb) :: translate      ! Flag to activate/deactivate translations
  logical(lb) :: rotate         ! Flag to activate/deactivate rotations
  logical(lb) :: computeProps   ! Flag to activate/deactivate energy computations
  integer(ib) :: rotationMode   ! Algorithm used for free rotation of rigid bodies
end type tOpts

type, bind(C) :: tEmDee
  integer(ib) :: builds         ! Number of neighbor list builds
  real(rb)    :: pairTime       ! Time taken in force calculations
  real(rb)    :: totalTime      ! Total time since initialization
  real(rb)    :: Potential      ! Total potential energy of the system
  real(rb)    :: Kinetic        ! Total kinetic energy of the system
  real(rb)    :: Rotational     ! Rotational kinetic energy of the system
  real(rb)    :: Virial         ! Total internal virial of the system
  integer(ib) :: DOF            ! Total number of degrees of freedom
  integer(ib) :: rotationDOF    ! Number of rotational degrees of freedom
  logical(lb) :: UpToDate       ! Flag to attest whether energies have been computed
  type(c_ptr) :: Data           ! Pointer to system data
  type(tOpts) :: Options        ! List of options to change EmDee's behavior
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
    type(pairModelContainer) :: noPair
    type(modelContainer)     :: noCoul

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
    me%threadAtoms = (N + threads - 1)/threads

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
    allocate( me%cell(0) )
    allocate( me%atomCell(N) )

    ! Allocate variables related to rigid bodies:
    call allocate_rigid_bodies( me, bodies )

    ! Allocate memory for list of atoms per cell:
    call me % cellAtom % allocate( N, 0 )

    ! Allocate memory for neighbor lists:
    allocate( me%neighbor(threads) )
    call me % neighbor % allocate( extra, N )

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
    allocate( me%layer_energy(me%nlayers), source = zero )
    allocate( me%layer_virial(me%nlayers), source = zero )

    ! Set up mutable entities:
    EmDee_system % builds = 0
    EmDee_system % pairTime = zero
    EmDee_system % totalTime = zero
    EmDee_system % Potential = zero
    EmDee_system % Kinetic = zero
    EmDee_system % Rotational = zero
    EmDee_system % DOF = 3*(N - 1)
    EmDee_system % rotationDOF = 0
    EmDee_system % data = c_loc(me)
    EmDee_system % Options % translate = .true.
    EmDee_system % Options % rotate = .true.
    EmDee_system % Options % computeProps = .true.
    EmDee_system % Options % rotationMode = 0

  end function EmDee_system

!===================================================================================================

  subroutine EmDee_switch_model_layer( md, layer ) bind(C,name="EmDee_switch_model_layer")
    type(tEmDee), value :: md
    integer(ib),  value :: layer

    type(tData), pointer :: me

    call c_f_pointer( md%data, me )
    if (layer /= me%layer) then
      if ((layer < 1).or.(layer > me%nlayers)) call error( "switch_model_layer", "out of range" )
      if (me%initialized) call update_forces( md, layer )
      me%layer = layer
    end if

  end subroutine EmDee_switch_model_layer

!===================================================================================================

  subroutine EmDee_set_pair_model( md, itype, jtype, model, kCoul ) &
    bind(C,name="EmDee_set_pair_model")
    type(tEmDee), value :: md
    integer(ib),  value :: itype, jtype
    type(c_ptr),  value :: model
    real(rb),     value :: kCoul

    character(*), parameter :: task = "pair model setting"

    integer :: layer, ktype
    type(tData), pointer :: me
    type(pairModelContainer), pointer :: container

    call c_f_pointer( md%data, me )
    if (.not.ranged( [itype,jtype], me%ntypes )) then
      call error( task, "provided type index is out of range" )
    end if
    if (me%initialized) call error( task, "cannot set model after coordinates have been defined" )
    if (.not.c_associated(model)) call error( task, "a valid pair model must be provided" )
    call c_f_pointer( model, container )
    select type ( pair_model => container%model )
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

    character(*), parameter :: task = "pair multimodel setting"

    integer :: layer, ktype
    character(5) :: C
    type(tData), pointer :: me
    type(pairModelContainer), pointer :: container

    call c_f_pointer( md%data, me )
    if (.not.ranged( [itype,jtype], me%ntypes )) then
      call error( task, "provided type index is out of range" )
    end if

    write(C,'(I5)') me%nlayers
    C = adjustl(C)

    if (me%initialized) call error( task, "cannot set model after coordinates have been defined" )
    do layer = 1, me%nlayers
      if (c_associated(model(layer))) then
        call c_f_pointer( model(layer), container )
      else
        call error( task, trim(C)//" valid pair models must be provided" )
      end if
      select type (pair_model => container%model )
        class is (cPairModel)
          call set_pair_type( me, itype, jtype, layer, container, kCoul(layer) )
        class default
          call error( task, trim(C)//" valid pair models must be provided" )
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

  subroutine EmDee_set_coul_model( md, model ) bind(C,name="EmDee_set_coul_model")
    type(tEmDee), value :: md
    type(c_ptr),  value :: model

    character(*), parameter :: task = "coulomb model setting"

    integer :: layer
    type(tData), pointer :: me
    type(modelContainer), pointer :: container

    call c_f_pointer( md%data, me )
    call c_f_pointer( model, container )

    if (me%initialized) call error( task, "cannot set model after coordinates have been defined" )
    if (.not.c_associated(model)) call error( task, "a valid coulomb model must be provided" )
    do layer = 1, me%nlayers
      me%coul(layer) = container
      select type ( coul_model => me%coul(layer)%model )
        class is (cCoulModel)
          call coul_model % shifting_setup( me%Rc )
        class default
          call error( task, "a valid coulomb model must be provided" )
      end select
    end do
    me % multilayer_coulomb = .false.

  end subroutine EmDee_set_coul_model

!===================================================================================================

  subroutine EmDee_set_coul_multimodel( md, model ) bind(C,name="EmDee_set_coul_multimodel")
    type(tEmDee), value :: md
    type(c_ptr),  intent(in) :: model(*)

    character(*), parameter :: task = "coulomb multimodel setting"

    integer :: layer
    character(5) :: C
    type(tData), pointer :: me
    type(modelContainer), pointer :: container

    call c_f_pointer( md%data, me )

    write(C,'(I5)') me%nlayers
    C = adjustl(C)

    if (me%initialized) call error( task, "cannot set model after coordinates have been defined" )
    do layer = 1, me%nlayers
      if (.not.c_associated(model(layer))) then
        call error( task, trim(C)//" valid coulomb models must be provided" )
      end if
      call c_f_pointer( model(layer), container )
      me%coul(layer) = container
      select type ( coul_model => me%coul(layer)%model )
        class is (cCoulModel)
          call coul_model % shifting_setup( me%Rc )
        class default
          call error( task, trim(C)//" valid coulomb models must be provided" )
      end select
    end do
    me % multilayer_coulomb = .true.

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

  subroutine EmDee_upload( md, option, address ) bind(C,name="EmDee_upload")
    type(tEmDee),      intent(inout) :: md
    character(c_char), intent(in)    :: option(*)
    type(c_ptr),       value         :: address

    real(rb) :: twoKEt, twoKEr
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
        me%invL = one/scalar
        me%invL2 = me%invL**2
        if (.not.me%initialized) then
          me%initialized = allocated( me%R )
          if (me%initialized) then
            call check_actual_interactions( me )
            !$omp parallel num_threads(me%nthreads)
            call update_rigid_bodies( me, omp_get_thread_num() + 1 )
            !$omp end parallel
            md%DOF = 3*me%nfree + sum(me%body%dof) - 3
          end if
        end if
        if (me%initialized) call compute_forces( md )

      case ("coordinates")
        if (.not.allocated( me%R )) allocate( me%R(3,me%natoms) )
        call c_f_pointer( address, Matrix, [3,me%natoms] )
        !$omp parallel num_threads(me%nthreads)
        call assign_coordinates( me, omp_get_thread_num() + 1, Matrix )
        !$omp end parallel
        if (.not.me%initialized) then
          me%initialized = me%Lbox > zero
          if (me%initialized) then
            call check_actual_interactions( me )
            !$omp parallel num_threads(me%nthreads)
            call update_rigid_bodies( me, omp_get_thread_num() + 1 )
            !$omp end parallel
            md%DOF = 3*me%nfree + sum(me%body%dof) - 3
          end if
        end if
        if (me%initialized) call compute_forces( md )

      case ("momenta")
        if (.not.me%initialized) call error( "upload", "box and coordinates have not been defined" )
        call c_f_pointer( address, Matrix, [3,me%natoms] )
        !$omp parallel num_threads(me%nthreads) reduction(+:TwoKEt,TwoKEr)
        call assign_momenta( me, omp_get_thread_num() + 1, Matrix, twoKEt, twoKEr )
        !$omp end parallel
        md%Kinetic = half*(twoKEt + twoKEr)
        md%Rotational = half*twoKEr

      case ("forces")
        if (.not.me%initialized) call error( "upload", "box and coordinates have not been defined" )
        call c_f_pointer( address, Matrix, [3,me%natoms] )
        !$omp parallel num_threads(me%nthreads)
        call assign_forces( omp_get_thread_num() + 1, Matrix )
        !$omp end parallel

      case ("charges")
        if (me%initialized) &
          call error( "upload", "cannot set charges after box and coordinates initialization" ) 
        call c_f_pointer( address, Vector, [me%natoms] )
        !$omp parallel num_threads(me%nthreads)
        call assign_charges( omp_get_thread_num() + 1, Vector )
        !$omp end parallel
        if (me%initialized) call compute_forces( md )

      case default
        call error( "upload", "invalid option" )

    end select

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine assign_forces( thread, Fext )
        integer,  intent(in) :: thread
        real(rb), intent(in) :: Fext(3,me%natoms)
        integer :: i, j
        do j = (thread - 1)*me%threadFreeAtoms + 1, min(thread*me%threadFreeAtoms, me%nfree)
          i = me%free(j)
          me%F(:,i) = Fext(:,i)
        end do
        do i = (thread - 1)*me%threadBodies + 1, min(thread*me%threadBodies, me%nbodies)
          call me % body(i) % force_torque_virial( Fext )
        end do
      end subroutine assign_forces
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine assign_charges( thread, Qext )
        integer,  intent(in) :: thread
        real(rb), intent(in) :: Qext(me%natoms)
        integer :: first, last
        first = (thread - 1)*me%threadAtoms + 1
        last = min(thread*me%threadAtoms, me%natoms)
        me%charge(first:last) = Qext(first:last)
        me%charged(first:last) = Qext(first:last) /= zero
      end subroutine assign_charges
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine assign_coordinates( me, thread, Rext )
        type(tData), intent(inout) :: me
        integer,     intent(in)    :: thread
        real(rb),    intent(in)    :: Rext(3,me%natoms)
        integer :: first, last
        first = (thread - 1)*me%threadAtoms + 1
        last = min(thread*me%threadAtoms, me%natoms)
        me%R(:,first:last) = Rext(:,first:last)
      end subroutine assign_coordinates
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
          end associate
        end do
      end subroutine update_rigid_bodies
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
        !$omp parallel num_threads(me%nthreads)
        call download( omp_get_thread_num() + 1, me%F, matrix )
        !$omp end parallel

      case ("multienergy")
        if (.not.md%UpToDate) call error( "download", "layer energies are outdated" )
        call c_f_pointer( address, vector, [me%nlayers] )
        vector = me%layer_energy

      case ("multivirial")
        if (.not.md%UpToDate) call error( "download", "layer virials are outdated" )
        call c_f_pointer( address, vector, [me%nlayers] )
        vector = me%layer_virial

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
      subroutine download( thread, origin, destiny )
        integer,  intent(in)    :: thread
        real(rb), intent(in)    :: origin(3,me%natoms)
        real(rb), intent(inout) :: destiny(3,me%natoms)
        integer :: first, last
        first = (thread - 1)*me%threadAtoms + 1
        last = min(thread*me%threadAtoms, me%natoms)
        destiny(:,first:last) = origin(:,first:last)
      end subroutine download
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
        if (.not.me%initialized) call error( "random_momenta", "coordinates have not defined" )
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

!  subroutine EmDee_save_state( md, rigid )
!    type(tEmDee), intent(inout) :: md
!    integer(ib),  intent(in)    :: rigid
!    if (rigid /= 0) then
!    else
!    end if
!  end subroutine EmDee_save_state

!===================================================================================================

!  subroutine EmDee_restore_state( md )
!    type(tEmDee), intent(inout) :: md
!  end subroutine EmDee_restore_state

!===================================================================================================

  subroutine EmDee_boost( md, lambda, alpha, dt ) bind(C,name="EmDee_boost")
    type(tEmDee), intent(inout) :: md
    real(rb),     value         :: lambda, alpha, dt

    real(rb) :: CP, CF, Ctau, twoKEt, twoKEr, KEt
    type(tData), pointer :: me

    call c_f_pointer( md%data, me )

    CF = phi(alpha*dt)*dt
    CP = one - alpha*CF
    CF = lambda*CF
    Ctau = two*CF

    twoKEt = zero
    twoKEr = zero
    !$omp parallel num_threads(me%nthreads) reduction(+:twoKEt,twoKEr)
    call boost( omp_get_thread_num() + 1, twoKEt, twoKEr )
    !$omp end parallel
    if (md%Options%translate) then
      KEt = half*twoKEt
    else
      KEt = md%Kinetic - md%Rotational
    end if
    if (md%Options%rotate) md%Rotational = half*twoKEr
    md%Kinetic = KEt + md%Rotational

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine boost( thread, twoKEt, twoKEr )
        integer,  intent(in)    :: thread
        real(rb), intent(inout) :: twoKEt, twoKEr
        integer  :: i, j
        do i = (thread - 1)*me%threadBodies + 1, min(thread*me%threadBodies, me%nbodies)
          associate(b => me%body(i))
            if (md%Options%translate) then
              b%pcm = CP*b%pcm + CF*b%F
              twoKEt = twoKEt + b%invMass*sum(b%pcm*b%pcm)
            end if
            if (md%Options%rotate) then
              call b%assign_momenta( CP*b%pi + matmul( matrix_C(b%q), Ctau*b%tau ) )
              twoKEr = twoKEr + sum(b%MoI*b%omega*b%Omega)
            end if
          end associate
        end do
        if (md%Options%translate) then
          do i = (thread - 1)*me%threadFreeAtoms + 1, min(thread*me%threadFreeAtoms, me%nfree)
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
            if (md%Options%translate) b%rcm = cR*b%rcm + cP*b%invMass*b%pcm
            if (md%Options%rotate) then
              if (md%Options%rotationMode == 0) then
               call b % rotate_exact( dt )
              else
                call b % rotate_no_squish( dt, n = md%options%rotationMode )
              end if
              forall (j=1:3) me%R(j,b%index) = b%rcm(j) + b%delta(j,:)
            end if
          end associate
        end do
        if (md%Options%translate) then
          do i = (thread - 1)*me%threadFreeAtoms + 1, min(thread*me%threadFreeAtoms, me%nfree)
            j = me%free(i)
            me%R(:,j) = cR*me%R(:,j) + cP*me%P(:,j)*me%invMass(j)
          end do
        end if
      end subroutine move
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine EmDee_move

!===================================================================================================
!                              A U X I L I A R Y   P R O C E D U R E S
!===================================================================================================

  subroutine compute_forces( md )
    type(tEmDee), intent(inout) :: md

    integer  :: M
    real(rb) :: time, E, W
    logical  :: buildList
    real(rb), allocatable :: Rs(:,:), Fs(:,:,:), Elayer(:,:), Wlayer(:,:)
    type(tData),  pointer :: me

    call c_f_pointer( md%data, me )
    md%pairTime = md%pairTime - omp_get_wtime()

    allocate( Rs(3,me%natoms), Fs(3,me%natoms,me%nthreads) )
    allocate( Elayer(me%nlayers,me%nthreads), Wlayer(me%nlayers,me%nthreads) )
    Rs = me%invL*me%R

    buildList = maximum_approach_sq( me%natoms, me%R - me%R0 ) > me%skinSq
    if (buildList) then
      M = floor(ndiv*me%Lbox/me%xRc)
      call distribute_atoms( me, max(M,2*ndiv+1), Rs )
      me%R0 = me%R
      md%builds = md%builds + 1
    endif

    md%UpToDate = md%Options%computeProps
    !$omp parallel num_threads(me%nthreads) reduction(+:E,W)
    block
      integer :: thread
      thread = omp_get_thread_num() + 1
      associate( F => Fs(:,:,thread), EL => Elayer(:,thread), WL => Wlayer(:,thread) )
        if (buildList) then
          call find_pairs_and_compute( me, thread, md%UpToDate, Rs, F, E, W, EL, WL )
        else
          call compute_pairs( me, thread, md%UpToDate, Rs, F, E, W, EL, WL )
        end if
        if (me%bonds%exist) call compute_bonds( me, thread, Rs, F, E, W )
        if (me%angles%exist) call compute_angles( me, thread, Rs, F, E, W )
        if (me%dihedrals%exist) call compute_dihedrals( me, thread, Rs, F, E, W )
      end associate
    end block
    !$omp end parallel

    me%F = me%Lbox*sum(Fs,3)
    md%Virial = W
    if (me%nbodies /= 0) call rigid_body_forces( me, md%Virial )

    if (md%UpToDate) then
      md%Potential = E
      me%layer_energy = sum(Elayer,2)
      me%layer_virial = sum(Wlayer,2)
    end if

    time = omp_get_wtime()
    md%pairTime = md%pairTime + time
    md%totalTime = time - me%startTime

  end subroutine compute_forces

!===================================================================================================

  subroutine update_forces( md, layer )
    type(tEmDee), intent(inout) :: md
    integer,      intent(in)    :: layer

    real(rb) :: DE, DW, time
    real(rb), allocatable :: Rs(:,:), DFs(:,:,:)
    type(tData),  pointer :: me

    call c_f_pointer( md%data, me )
    md%pairTime = md%pairTime - omp_get_wtime()

    allocate( Rs(3,me%natoms), DFs(3,me%natoms,me%nthreads) )
    Rs = me%invL*me%R

    !$omp parallel num_threads(me%nthreads) reduction(+:DE,DW)
    block
      integer :: thread
      thread = omp_get_thread_num() + 1
      call update_pairs( me, thread, Rs, DFs(:,:,thread), DE, DW, layer )
    end block
    !$omp end parallel

    me%F = me%F + me%Lbox*sum(DFs,3)
    md%Potential = md%Potential + DE
    md%Virial = md%Virial + DW
    if (me%nbodies /= 0) call rigid_body_forces( me, md%Virial )

    time = omp_get_wtime()
    md%pairTime = md%pairTime + time
    md%totalTime = time - me%startTime

  end subroutine update_forces

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

  character(sl) function string( carray )
    character(c_char), intent(in) :: carray(*)
    integer :: i
    string = ""
    do i = 1, sl
      if (carray(i) == c_null_char) return
      string(i:i) = carray(i)
    end do
  end function string

!===================================================================================================

end module EmDeeCode
