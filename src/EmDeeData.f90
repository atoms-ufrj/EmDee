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

type, private :: tCell
  integer :: neighbor(nbcells)
end type tCell

type, bind(C), public :: tEnergy
  real(rb)    :: Potential            ! Total potential energy of the system
  real(rb)    :: Dispersion           ! Dispersion (vdW) part of the potential energy
  real(rb)    :: Coulomb              ! Electrostatic part of the potential energy
  real(rb)    :: Bond
  real(rb)    :: Angle
  real(rb)    :: Dihedral
  real(rb)    :: ShadowPotential
  logical(lb) :: UpToDate             ! Flag to attest whether energies have been computed
end type tEnergy

type, bind(C), public :: tKinetic
  real(rb)    :: Total                ! Total kinetic energy of the system
  real(rb)    :: TransPart(3)         ! Translational kinetic energy at each dimension
  real(rb)    :: Rotational           ! Total rotational kinetic energy of the system
  real(rb)    :: RotPart(3)           ! Rotational kinetic energy around each principal axis
  real(rb)    :: ShadowKinetic
  real(rb)    :: ShadowRotational
  logical(lb) :: UpToDate             ! Flag to attest whether energies have been computed
end type tKinetic

type, bind(C), public :: tVirial
  real(rb)    :: Total                ! Total internal virial of the system
  real(rb)    :: Body                 ! Rigid body contribution to the internal virial
end type tVirial

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

  real(rb), pointer :: Lbox => null()

  real(rb) :: Rc                          ! Cut-off distance
  real(rb) :: skin                        ! Neighbor list skin width
  real(rb) :: RcSq                        ! Cut-off distance squared
  real(rb) :: xRc                         ! Extended cutoff distance (including skin)
  real(rb) :: xRcSq                       ! Extended cutoff distance squared
  real(rb) :: skinSq                      ! Square of the neighbor list skin width
  real(rb) :: InRc                        ! Internal cut-off distance
  real(rb) :: InRcSq                      ! Internal cut-off distance squared
  real(rb) :: xInRcSq                     ! Extended internal cutoff distance squared
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
  integer, allocatable :: atomBody(:)     ! Array containing the rigid body that contains each atom
  integer, allocatable :: free(:)         ! Pointer to the list of independent atoms
  integer, allocatable :: atomsInCell(:)  ! Number of atoms in each cell

  real(rb),    pointer :: R(:,:) => null()   ! Coordinates of all atoms
  real(rb),    pointer :: P(:,:) => null()   ! Momenta of all atoms
  type(tBody), pointer :: body(:) => null()  ! Pointer to the rigid bodies present in the system

  real(rb), pointer     :: F(:,:)         ! Resultant forces on all atoms
  real(rb), allocatable :: charge(:)      ! Electric charges of all atoms
  real(rb), allocatable :: mass(:)        ! Masses of all atoms
  real(rb), allocatable :: invMass(:)     ! Inverses of atoms masses
  real(rb), allocatable :: R0(:,:)        ! Position of each atom at latest neighbor list building
  logical,  allocatable :: charged(:)     ! Flag to determine if a particle is charged

  type(tCell), allocatable :: cell(:)      ! Array containing all neighbor cells of each cell
  type(tList), allocatable :: neighbor(:)  ! Pointer to neighbor lists

  type(pairContainer), allocatable :: pair(:,:,:)
  type(coulContainer), allocatable :: coul(:)
  class(cKspaceModel), allocatable :: kspace

  logical,  allocatable :: multilayer(:,:)
  logical,  allocatable :: overridable(:,:)
  logical,  allocatable :: interact(:,:)
  logical,  allocatable :: pairs_exist(:)

  ! Fields related to model layers:
  logical,       allocatable :: bonded(:)
  logical,       allocatable :: useInRc(:)
  logical,       allocatable :: forcesUpToDate(:)
  real(rb),      allocatable :: layerF(:,:,:)
  type(tEnergy), allocatable :: layerEnergy(:)
  type(tVirial), allocatable :: layerVirial(:)

  logical :: multilayer_coulomb
  logical :: kspace_active

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
    logical :: interact(me%ntypes,me%ntypes,me%nlayers)

    interact = .false.
    do i = 1, me%ntypes
      neutral(i) = count((me%atomType == i) .and. me%charged) == 0
      do j = 1, i
        associate( pair => me%pair(i,j,:) )
          do k = 1, me%nlayers
            no_pair = same_type_as( pair(k)%model, PairNone )
            no_coul = same_type_as( me%coul(k)%model, CoulNone ) .or. (.not.pair(k)%coulomb)
            coul_only = no_pair .and. (.not.no_coul)
            inert = (no_pair .and. no_coul) .or. (coul_only .and. neutral(i) .and. neutral(j))
            interact(i,j,k) = .not.inert
          end do
          me%interact(i,j) = any(interact(i,j,:))
          me%interact(j,i) = me%interact(i,j)
        end associate
      end do
    end do
    do k = 1, me%nlayers
      me%pairs_exist(k) = any(interact(:,:,k))
    end do

  end subroutine check_actual_interactions

!===================================================================================================

  subroutine set_pair_type( me, itype, jtype, layer, container, kCoul )
    type(tData),          intent(inout) :: me
    integer(ib),          intent(in)    :: itype, jtype, layer
    type(modelContainer), intent(in)    :: container
    real(rb),             intent(in)    :: kCoul

    integer :: ktype
    real(rb) :: layerRc(me%nlayers)

    layerRc = merge(me%InRc, me%Rc, me%useInRc)
    select type (pmodel => container%model)
      class is (cPairModel)
        associate (pair => me%pair(:,:,layer))
          if (itype == jtype) then
            pair(itype,itype) = container
            pair(itype,itype)%coulomb = kCoul /= zero
            if (pair(itype,itype)%coulomb) pair(itype,itype)%kCoul = kCoul
            call pair(itype,itype) % model % modifier_setup( layerRc(layer) )
            do ktype = 1, me%ntypes
              if ((ktype /= itype).and.me%overridable(itype,ktype)) then
                pair(itype,ktype) = pair(ktype,ktype) % mix( pair(itype,itype) )
                call pair(itype,ktype) % model % modifier_setup( layerRc(layer) )
                pair(ktype,itype) = pair(itype,ktype)
              end if
            end do
          else
            pair(itype,jtype) = container
            pair(itype,jtype)%coulomb = kCoul /= zero
            if (pair(itype,jtype)%coulomb) pair(itype,jtype)%kCoul = kCoul
            call pair(itype,jtype) % model % modifier_setup( layerRc(layer) )
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
        ! Otherwise, 1 <= index(i) <= number of distinct positive entries which are non-unique

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

    integer :: i, bodyDoF
    real(rb) :: L(3), kspaceRc, layerRc(me%nlayers)
    integer, allocatable :: bodyAtom(:), coulomb_with_kspace(:)
    logical, allocatable :: kspace_required(:)

    !$omp parallel num_threads(me%nthreads)
    call update_rigid_bodies( me, omp_get_thread_num() + 1 )
    !$omp end parallel

    bodyDoF = sum(me%body%dof)
    RotDoF = bodyDoF - 3*me%nbodies
    DoF = 3*me%nfree + bodyDoF - 3

    call check_actual_interactions( me )

    ! Find out which coulomb models require a kspace solver:
    kspace_required = [(me%coul(i)%model%requires_kspace, i=1, me%nlayers)]
    coulomb_with_kspace = pack([(i, i=1, me%nlayers)], kspace_required)

    ! Check if all coulomb models requiring kspace have the same cutoff:
    if (size(coulomb_with_kspace) > 0) then
      layerRc = merge(me%InRc, me%Rc, me%useInRc)
      kspaceRc = layerRc(coulomb_with_kspace(1))
      if (any(layerRc(coulomb_with_kspace) /= kspaceRc)) then
        call error(task, "all layers with ewald-like coulomb models must have the same cutoff")
      end if
    end if

    ! Now check if demand and supply are inconsistent:
    if (any(kspace_required) .neqv. me%kspace_active) then
      if (me%kspace_active) then
        ! Deactivate kspace solver since it is not required:
        me%kspace_active = .false.
      else
        call error( task, "a kspace solver is required, but has not been defined" )
      end if
    end if

    ! Initialize kspace solver and related coulomb models, if needed:
    if (me%kspace_active) then
      bodyAtom = [(me%body(i)%index, i=1, me%nbodies)]
      L = me%Lbox
      call me % kspace % initialize( me%nthreads, kspaceRc, L, me%atomType, me%charge, &
                                     me%R, me%body%NP, bodyAtom, me%pair%kCoul )
      do i = 1, size(coulomb_with_kspace)
        call me % coul(coulomb_with_kspace(i)) % model % kspace_setup( me%kspace%alpha )
      end do
    end if

    me%initialized = .true.

  end subroutine perform_initialization

!===================================================================================================

  subroutine update_rigid_bodies( me, thread )
    type(tData),  intent(inout) :: me
    integer,      intent(in)    :: thread

    integer :: i, j
    real(rb) :: L, invL
    real(rb), allocatable :: R(:,:)

    L = me%Lbox
    invL = one/L
    do j = (thread - 1)*me%threadBodies + 1, min(thread*me%threadBodies, me%nbodies)
      associate(b => me%body(j))
        R = me%R(:,b%index)
        forall (i = 2:b%NP) R(:,i) = R(:,i) - L*anint(invL*(R(:,i) - R(:,1)))
        call b % update( R )
        me%R(:,b%index) = R
      end associate
    end do

  end subroutine update_rigid_bodies

!===================================================================================================

  subroutine compute_bonds( me, threadId, R, F, Potential, Virial, Ecoul )
    type(tData), intent(inout) :: me
    integer,     intent(in)    :: threadId
    real(rb),    intent(in)    :: R(3,me%natoms)
    real(rb),    intent(inout) :: F(3,me%natoms), Potential, Virial, Ecoul

    integer  :: m, nbonds, itype, jtype
    real(rb) :: invL2, invR2, E, W, EL, WL, QiQj
    real(rb) :: Rij(3), Fij(3)
    logical  :: bonded, kspace

    bonded = me%bonded(me%layer)
    kspace = me%coul(me%layer)%model%requires_kspace
    if (me%bonds%exist.and.(bonded.or.kspace)) then
      nbonds = (me%bonds%number + me%nthreads - 1)/me%nthreads
      invL2 = one/me%Lbox**2
      do m = (threadId - 1)*nbonds + 1, min( me%bonds%number, threadId*nbonds )
        associate(bond => me%bonds%item(m))
          Rij = R(:,bond%i) - R(:,bond%j)
          Rij = Rij - anint(Rij)
          invR2 = invL2/sum(Rij*Rij)
          if (bonded) then
            select type (model => bond%model)
              include "compute_bond.f90"
            end select
            Potential = Potential + E
            Virial = Virial + W
          else
            W = zero
          end if
          if (kspace) then
            itype = me%atomType(bond%i)
            jtype = me%atomType(bond%j)
            QiQj = me%pair(itype,jtype,me%layer)%kCoul*me%charge(bond%i)*me%charge(bond%j)
            call me % kspace % discount( EL, WL, one/invR2, QiQj )
            Ecoul = Ecoul + EL
            W = W + WL
          end if
          Fij = W*invR2*me%Lbox*Rij
          F(:,bond%i) = F(:,bond%i) + Fij
          F(:,bond%j) = F(:,bond%j) - Fij
        end associate
      end do
    end if
  end subroutine compute_bonds

!===================================================================================================

  subroutine compute_angles( me, threadId, R, F, Potential, Virial, Ecoul )
    type(tData), intent(inout) :: me
    integer,     intent(in)    :: threadId
    real(rb),    intent(in)    :: R(3,me%natoms)
    real(rb),    intent(inout) :: F(3,me%natoms), Potential, Virial, Ecoul

    integer  :: i, j, k, m, nangles, itype, ktype
    real(rb) :: aa, bb, ab, theta, Ea, Fa, factor, QiQk
    real(rb) :: Fi(3), Fk(3), avec(3), bvec(3)
    logical  :: bonded, kspace

    bonded = me%bonded(me%layer)
    kspace = me%coul(me%layer)%model%requires_kspace
    if (me%angles%exist.and.(bonded.or.kspace)) then
      nangles = (me%angles%number + me%nthreads - 1)/me%nthreads
      do m = (threadId - 1)*nangles + 1, min( me%angles%number, threadId*nangles )
        associate (angle => me%angles%item(m))
          i = angle%i
          j = angle%j
          k = angle%k
          avec = R(:,i) - R(:,j)
          bvec = R(:,k) - R(:,j)
          avec = me%Lbox*(avec - anint(avec))
          bvec = me%Lbox*(bvec - anint(bvec))
          if (bonded) then
            aa = sum(avec*avec)
            bb = sum(bvec*bvec)
            ab = sum(avec*bvec)
            theta = acos(ab/sqrt(aa*bb))
            select type (model => angle%model)
              include "compute_angle.f90"
            end select
            factor = Fa/sqrt(aa*bb - ab*ab)
            Fi = ((ab/aa)*avec - bvec)*factor
            Fk = ((ab/bb)*bvec - avec)*factor
            F(:,i) = F(:,i) + Fi
            F(:,k) = F(:,k) + Fk
            F(:,j) = F(:,j) - (Fi + Fk)
            Potential = Potential + Ea
            Virial = Virial + sum(Fi*avec + Fk*bvec)
          end if
          if (kspace) then
            block
              real(rb) :: Rik(3), RikSq, Fik(3), EL, WL
              Rik = avec - bvec
              RikSq = sum(Rik*Rik)
              itype = me%atomType(i)
              ktype = me%atomType(k)
              QiQk = me%pair(itype,ktype,me%layer)%kCoul*me%charge(i)*me%charge(k)
              call me % kspace % discount( EL, WL, RikSq, QiQk )
              Ecoul = Ecoul + EL
              Fik = WL*Rik/RikSq
              F(:,i) = F(:,i) + Fik
              F(:,k) = F(:,k) - Fik
            end block
          end if
        end associate
      end do
    end if
  end subroutine compute_angles

!===================================================================================================

  subroutine compute_dihedrals( me, threadId, R, F, Potential, Virial )
    type(tData), intent(inout) :: me
    integer,     intent(in)    :: threadId
    real(rb),    intent(in)    :: R(3,me%natoms)
    real(rb),    intent(inout) :: F(3,me%natoms), Potential, Virial

    integer  :: m, ndihedrals, i, j
    real(rb) :: Rc2, Ed, Fd, r2, invL2, invR2, invR, Eij, Wij, Qi, Qj
    real(rb) :: Rj(3), Rk(3), Fi(3), Fk(3), Fl(3), Fij(3)
    real(rb) :: normRkj, normX, a, b, phi, factor14
    real(rb) :: rij(3), rkj(3), rlk(3), x(3), y(3), z(3), u(3), v(3), w(3)

    if (.not.me%dihedrals%exist) return
    invL2 = one/me%Lbox**2
    Rc2 = merge(me%InRcSq, me%RcSq, me%useInRc(me%layer))*invL2
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
          if (r2 < Rc2) then
            invR2 = invL2/r2
            invR = sqrt(invR2)
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

  subroutine compute_pairs( me, thread, compute, Rs, F, Epair, Ecoul, Wpair, Wcoul )
    type(tData), intent(inout) :: me
    integer,     intent(in)    :: thread
    logical(lb), intent(in)    :: compute
    real(rb),    intent(in)    :: Rs(3,me%natoms)
    real(rb),    intent(out)   :: F(3,me%natoms), Epair, Ecoul, Wpair, Wcoul

    real(rb) :: Rc2, L2, invL2
    integer, pointer :: upper(:)

    F = zero
    Wpair = zero
    Wcoul = zero
    Epair = zero
    if (me%pairs_exist(me%layer)) then
      L2 = me%Lbox**2
      invL2 = one/L2
      if (me%useInRc(me%layer)) then
        Rc2 = me%InRcSq*invL2
        upper => me%neighbor(thread)%middle(:)
      else
        Rc2 = me%RcSq*invL2
        upper => me%neighbor(thread)%last(:)
      end if
      if (compute) then
#       define compute
#       include "compute.f90"
#       undef compute
      else
#       include "compute.f90"
      end if
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

  subroutine compute_kspace( me, Rs, Elong, F )
    type(tData), intent(inout) :: me
    real(rb),    intent(in)    :: Rs(3,me%natoms)
    real(rb),    intent(inout) :: Elong, F(3,me%natoms)

    if (me%coul(me%layer)%model%requires_kspace) then
      call me % kspace % prepare( Rs )
      call me % kspace % compute( me%layer, Elong, F )
      call me % kspace % discount_rigid_pairs( me%layer, me%Lbox*[1,1,1], me%R, Elong )
    end if

  end subroutine compute_kspace

!===================================================================================================
  !
  ! subroutine update_pairs( me, thread, Rs, DF, DEpair, DEcoul, DWpair, DWcoul, new_layer )
  !   type(tData), intent(in)  :: me
  !   integer,     intent(in)  :: thread, new_layer
  !   real(rb),    intent(in)  :: Rs(3,me%natoms)
  !   real(rb),    intent(out) :: DF(3,me%natoms), DEpair, DEcoul, DWpair, DWcoul
  !
  !   integer  :: i, j, k, m, itype, jtype, firstAtom, lastAtom
  !   real(rb) :: invL2, Rc2, r2, invR2, invR, Wij, Qi, QiQj, Wij, rFc, DWij, DWij
  !   real(rb) :: Rij(3), Ri(3), DFi(3), DFij(3)
  !   logical  :: icharged, ijcharged
  !   logical  :: participant(me%ntypes)
  !   integer,  allocatable :: jlist(:)
  !   real(rb), allocatable :: Rvec(:,:)
  !
  !   participant = any(me%multilayer,dim=2)
  !   invL2 = one/me%Lbox**2
  !   Rc2 = me%RcSq*invL2
  !
  !   DF = zero
  !   DWpair = zero
  !   DWcoul = zero
  !
  !   associate ( neighbor => me%neighbor(thread) )
  !     firstAtom = me%cellAtom%first(me%threadCell%first(thread))
  !     lastAtom = me%cellAtom%last(me%threadCell%last(thread))
  !     do k = firstAtom, lastAtom
  !       i = me%cellAtom%item(k)
  !       itype = me%atomType(i)
  !       if (participant(itype)) then
  !         Qi = me%charge(i)
  !         icharged = me%charged(i)
  !         Ri = Rs(:,i)
  !         DFi = zero
  !         jlist = neighbor%item(neighbor%first(i):neighbor%last(i))
  !         jlist = pack(jlist,me%multilayer(me%atomType(jlist),itype))
  !         Rvec = Rs(:,jlist)
  !         forall (m=1:size(jlist)) Rvec(:,m) = pbc(Ri - Rvec(:,m))
  !         do m = 1, size(jlist)
  !           j = jlist(m)
  !           Rij = Rvec(:,m)
  !           r2 = sum(Rij*Rij)
  !           if (r2 < Rc2) then
  !             invR2 = invL2/r2
  !             invR = sqrt(invR2)
  !             jtype = me%atomType(j)
  !             ijcharged = icharged.and.me%charged(j)
  !
  !             associate ( pair => me%pair(jtype,itype,new_layer) )
  !               associate ( model => pair%model )
  !                 select type ( model )
  !                   include "virial_compute_pair.f90"
  !                 end select
  !                 select case (model%modifier)
  !                   case (2) ! SHIFTED_FORCE
  !                     rFc = model%fshift/invR
  !                     Wij = Wij - rFc
  !                 end select
  !               end associate
  !               DWij = Wij
  !               if (ijcharged.and.pair%coulomb) then
  !                 QiQj = pair%kCoul*Qi*me%charge(j)
  !                 select type ( model => me%coul(me%layer)%model )
  !                   include "virial_compute_coul.f90"
  !                 end select
  !                 DWij = Wij
  !               else
  !                 DWij = zero
  !               end if
  !             end associate
  !
  !             associate ( pair => me%pair(jtype,itype,me%layer) )
  !               associate ( model => pair%model )
  !                 select type ( model )
  !                   include "virial_compute_pair.f90"
  !                 end select
  !                 select case (model%modifier)
  !                   case (2) ! SHIFTED_FORCE
  !                     rFc = model%fshift/invR
  !                     Wij = Wij - rFc
  !                 end select
  !               end associate
  !               DWij = DWij - Wij
  !               if (ijcharged.and.pair%coulomb) then
  !                 QiQj = pair%kCoul*Qi*me%charge(j)
  !                 select type ( model => me%coul(me%layer)%model )
  !                   include "virial_compute_coul.f90"
  !                 end select
  !                 DWij = DWij - Wij
  !               end if
  !             end associate
  !
  !             DFij = (DWij + DWij)*invR2*Rij
  !             DFi = DFi + DFij
  !             DF(:,j) = DF(:,j) - DFij
  !
  !             DWpair = DWpair + DWij
  !             DWcoul = DWcoul + DWij
  !
  !           end if
  !         end do
  !       end if
  !     end do
  !   end associate
  !   where (DF /= zero) DF = me%Lbox*DF
  !   DEpair = me%threadEpair(new_layer,thread) - me%threadEpair(me%layer,thread)
  !   DEcoul = me%threadEcoul(new_layer,thread) - me%threadEcoul(me%layer,thread)
  !
  !   contains
  !     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !     elemental function pbc( x )
  !       real(rb), intent(in) :: x
  !       real(rb)              :: pbc
  !       pbc = x - anint(x)
  !     end function pbc
  !     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! end subroutine update_pairs
  !
!===================================================================================================

  subroutine move( me, thread, R_factor, P_factor, dt, translate, rotate, mode )
    type(tData), intent(inout) :: me
    integer,     intent(in)    :: thread
    real(rb),    intent(in)    :: R_factor, P_factor, dt
    logical(lb), intent(in)    :: translate, rotate
    integer(ib), intent(in)    :: mode

    integer :: i, j, b1, bN

    b1 = (thread - 1)*me%threadBodies + 1
    bN = min(thread*me%threadBodies, me%nbodies)

    if (translate) then
      do i = b1, bN
        associate ( b => me%body(i) )
          b%rcm = R_factor*b%rcm + P_factor*b%invMass*b%pcm
        end associate
      end do
      do i = (thread - 1)*me%threadFreeAtoms + 1, min(thread*me%threadFreeAtoms, me%nfree)
        j = me%free(i)
        me%R(:,j) = R_factor*me%R(:,j) + P_factor*me%P(:,j)*me%invMass(j)
      end do
    end if

    if (rotate) then
      do i = b1, bN
        associate(b => me%body(i))
          if (mode == 0) then
           call b % rotate_exact( dt )
          else
            call b % rotate_no_squish( dt, n = mode )
          end if
          forall (j=1:3) me%R(j,b%index) = b%rcm(j) + b%delta(j,:)
        end associate
      end do
    end if

  end subroutine move

!===================================================================================================

  subroutine boost( me, thread, P_factor, F_factor, F, translate, rotate )
    type(tData), intent(inout) :: me
    integer,     intent(in)    :: thread
    real(rb),    intent(in)    :: P_factor, F_factor, F(3,me%natoms)
    logical(lb), intent(in)    :: translate, rotate

    integer  :: i, j, b1, bN
    real(rb) :: Ctau

    b1 = (thread - 1)*me%threadBodies + 1
    bN = min(thread*me%threadBodies, me%nbodies)
    do i = b1, bN
      call me % body(i) % force_and_torque( F )
    end do

    if (translate) then
      forall (i=b1:bN) me%body(i)%pcm = P_factor*me%body(i)%pcm + F_factor*me%body(i)%F
      do i = (thread - 1)*me%threadFreeAtoms + 1, min(thread*me%threadFreeAtoms, me%nfree)
        j = me%free(i)
        me%P(:,j) = P_factor*me%P(:,j) + F_factor*F(:,j)
      end do
    end if

    if (rotate) then
      Ctau = two*F_factor
      do i = b1, bN
        associate(b => me%body(i))
          call b % assign_momenta( P_factor*b%pi + matmul( matrix_C(b%q), Ctau*b%tau ) )
        end associate
      end do
    end if

  end subroutine boost

!===================================================================================================

  subroutine kinetic_energies( me, thread, translate, rotate, twoKEt, twoKEr )
    type(tData), intent(in)  :: me
    integer,     intent(in)  :: thread
    logical(lb), intent(in)  :: translate, rotate
    real(rb),    intent(out) :: twoKEt(3), twoKEr(3)
    integer :: i, j
    if (translate) then
      twoKEt = zero
      do i = (thread - 1)*me%threadBodies + 1, min(thread*me%threadBodies, me%nbodies)
        twoKEt = twoKEt + me%body(i)%invMass*me%body(i)%pcm**2
      end do
      do i = (thread - 1)*me%threadFreeAtoms + 1, min(thread*me%threadFreeAtoms, me%nfree)
        j = me%free(i)
        twoKEt = twoKEt + me%invMass(j)*me%P(:,j)**2
      end do
    end if
    if (rotate) then
      twoKEr = zero
      do i = (thread - 1)*me%threadBodies + 1, min(thread*me%threadBodies, me%nbodies)
        twoKEr = twoKEr + me%body(i)%MoI*me%body(i)%omega**2
      end do
    end if
  end subroutine kinetic_energies

!===================================================================================================

  function rigid_body_virial( me ) result( Virial )
    type(tData), intent(in) :: me
    real(rb)                :: Virial

    real(rb) :: W(me%nthreads)

    !$omp parallel num_threads(me%nthreads)
    block
      integer :: thread
      thread = omp_get_thread_num() + 1
      call partial_body_virial( thread, W(thread) )
    end block
    !$omp end parallel
    Virial = -sum(W)

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine partial_body_virial( thread, W )
        integer,  intent(in)  :: thread
        real(rb), intent(out) :: W
        integer :: i
        W = zero
        do i = (thread - 1)*me%threadBodies + 1, min(thread*me%threadBodies, me%nbodies)
          W = W + sum(me%F(:,me%body(i)%index)*me%body(i)%delta)
        end do
      end subroutine partial_body_virial
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end function rigid_body_virial

!===================================================================================================

end module EmDeeData
