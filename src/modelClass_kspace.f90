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

module kspaceModelClass

use global
use lists
use modelClass
use omp_lib
use math

implicit none

type tRigidPair
  integer(ib) :: i
  integer(ib) :: j
  integer(ib) :: itype
  integer(ib) :: jtype
  real(rb)    :: FbyR    ! Unscaled force-norm divided by distance
end type tRigidPair

!> An abstract class for kspace interaction models:
type, abstract, extends(cModel) :: cKspaceModel

  real(rb) :: alpha                               ! Ewald damping parameter alpha
  real(rb) :: beta                                ! Auxiliary constant: 2*alpha/sqrt(Pi)

  integer :: ntypes                               ! Number of distinc types of charged atoms
  integer,  allocatable :: type(:)                ! Distinct types of charged atoms

  type(tList) :: charge                           ! List of charged atoms and their charges
  integer,  allocatable :: ncharges(:)            ! Number of charged atoms of each type

  integer :: nTypePairs                           ! Number of types of atoms pairs
  real(rb), allocatable :: Erigid(:)              ! Unscaled energy of intrabody atom pair

  integer :: nRigidPairs                          ! Number of intrabody pairs (fixed distance)
  type(tRigidPair), allocatable :: rigidPair(:)   ! List of intrabody atom pairs

  integer :: nlayers                              ! Number of model layers
  real(rb), allocatable :: lambda(:,:,:)          ! Coulomb constants of type pairs
  real(rb), allocatable :: lambda1D(:,:)          ! Coulomb constants of type pairs (1D version)

  integer :: nthreads                             ! Number of parallel threads
  integer :: threadCharges                        ! Number of charges per thread
  integer :: threadRigidPairs                     ! Number of rigid pairs per thread
  integer :: threadTypePairs                      ! Number of type pairs per thread

  contains
    procedure :: initialize => cKspaceModel_initialize
    procedure :: discount_rigid_pairs => cKspaceModel_discount_rigid_pairs
    procedure :: discount => cKspaceModel_discount
    procedure(cKspaceModel_set_parameters), deferred :: set_parameters
    procedure(cKspaceModel_update),  deferred :: update
    procedure(cKspaceModel_prepare), deferred :: prepare
    procedure(cKspaceModel_compute), deferred :: compute
end type cKspaceModel

abstract interface

  subroutine cKspaceModel_set_parameters( me, Rc, L, alpha )
    import
    class(cKspaceModel), intent(inout) :: me
    real(rb),            intent(in)    :: Rc, L(3)
    real(rb),            intent(out)   :: alpha
  end subroutine cKspaceModel_set_parameters

  ! This procedure must be invoked whenever the box geometry and/or atom charges change:
  subroutine cKspaceModel_update( me, L )
    import
    class(cKspaceModel), intent(inout) :: me
    real(rb),            intent(in)    :: L(3)
  end subroutine cKspaceModel_update

  subroutine cKspaceModel_prepare( me, Rs )
    import
    class(cKspaceModel), intent(inout) :: me
    real(rb),            intent(in)    :: Rs(:,:)
  end subroutine cKspaceModel_prepare

  subroutine cKspaceModel_compute( me, layer, E, F )
    import
    class(cKspaceModel), intent(inout) :: me
    integer(ib),         intent(in)    :: layer
    real(rb),            intent(inout) :: E
    real(rb), optional,  intent(inout) :: F(:,:)
  end subroutine cKspaceModel_compute

end interface

contains

!---------------------------------------------------------------------------------------------------

  subroutine cKspaceModel_initialize( me, nthreads, Rc, L, types, charge, R, &
                                      groupSize, groupAtom, kCoul )
    class(cKspaceModel), intent(inout) :: me
    integer,             intent(in)    :: nthreads
    real(rb),            intent(in)    :: Rc, L(3)
    integer,             intent(in)    :: types(:)
    real(rb),            intent(in)    :: charge(size(types))
    integer,             intent(in)    :: groupSize(:), groupAtom(:)
    real(rb),            intent(in)    :: R(3,size(types)), kCoul(:,:,:)

    character(*), parameter :: task = "kspace model initialization"

    integer :: i, j, itype, jtype, ncharged, nChargedTypes
    logical :: charged(size(types))

    ! Mark and count charged atoms:
    charged = abs(charge) > epsilon(one)
    ncharged = count(charged)
    if (ncharged == 0) call error( task, "system has no charged atoms" )

    ! List distinct types of charged atoms:
    block
      integer, allocatable :: distinct(:), chargedTypes(:)

      chargedTypes = pack( types, charged )
      nChargedTypes = maxval(chargedTypes)
      allocate( distinct(nChargedTypes) )
      me%ntypes = 0
      do i = 1, size(chargedTypes)
        itype = chargedTypes(i)
        if (.not.any(distinct(1:me%ntypes) == itype)) then
          me%ntypes = me%ntypes + 1
          distinct(me%ntypes) = itype
        end if
      end do
      me%type = sorted(distinct(1:me%ntypes))
    end block

    ! Allocate and populate charge list:
    block
      integer :: first, last
      integer, allocatable :: seq(:), chargedAtoms(:)

      call me % charge % allocate( ncharged, me%ntypes, value = .true. )
      seq = [(i,i=1,size(types))]
      last = 0
      do i = 1, me%ntypes
        first = last + 1
        chargedAtoms = pack( seq, types == me%type(i) )
        last = last + size(chargedAtoms)
        me%charge%first(i) = first
        me%charge%last(i) = last
        me%charge%item(first:last) = chargedAtoms
        me%charge%value(first:last) = charge(chargedAtoms)
      end do
      me%ncharges = me%charge%last - me%charge%first + 1
    end block

    ! Set kspace model parameters:
    call me % set_parameters( Rc, L, me%alpha )
    me%beta = two*me%alpha/sqrt(Pi)

    ! List rigid pairs and precompute their contributions::
    block
      integer  :: maxnpairs, npairs, ii, i, jj, j, kk, k, prev, n
      integer  :: itype, jtype
      real(rb) :: invL(3), Ri(3), Rij(3), Qi, rsq, Eij, Wij

      integer,          allocatable :: local(:)
      type(tRigidPair), allocatable :: pair(:)

      me%nTypePairs = me%ntypes*(me%ntypes - 1)/2 + me%ntypes
      allocate( me%Erigid(me%nTypePairs), source = zero )

      maxnpairs = sum(groupSize*(groupSize-1)/2)
      allocate( pair(maxnpairs) )

      allocate( local(nChargedTypes) )
      local(me%type) = [(i,i=1,me%ntypes)]

      npairs = 0
      prev = 0
      invL = one/L
      do kk = 1, size(groupSize)
        n = groupSize(kk)
        do ii = 1, n - 1
          i = groupAtom(prev+ii)
          if (charged(i)) then
            Ri = R(:,i)
            Qi = charge(i)
            itype = local(types(i))
            do jj = ii + 1, n
              j = groupAtom(prev+jj)
              if (charged(j)) then
                jtype = local(types(j))

                Rij = Ri - R(:,j)
                Rij = Rij - L*anint(invL*Rij)
                rsq = sum(Rij**2)
                call me % discount( Eij, Wij, rsq, Qi*charge(j) )

                k = symm1D( itype, jtype )
                me%Erigid(k) = me%Erigid(k) + Eij

                npairs = npairs + 1
                pair(npairs) = tRigidPair( i, j, itype, jtype, Wij/rsq )
              end if
            end do
          end if
        end do
        prev = prev + n
      end do

      ! Discount self-energy:
      do itype = 1, me%ntypes
        k = symm1D( itype, itype )
        associate( Q => me%charge%value(me%charge%first(itype):me%charge%last(itype)) )
          me%Erigid(k) = me%Erigid(k) - me%alpha*sum(Q**2)/sqrt(Pi)
        end associate
      end do

      me%nRigidPairs = npairs
      me%rigidPair = pair(1:npairs)
    end block

    ! Store coulomb constants in both 1D and 2D formats:
    me%nlayers = size(kCoul,3)
    allocate( me%lambda(me%ntypes,me%ntypes,me%nlayers) )
    allocate( me%lambda1D(me%nTypePairs,me%nlayers) )
    do i = 1, me%ntypes
      itype = me%type(i)
      do j = i, me%ntypes
        jtype = me%type(j)
        associate ( lambda => kCoul(itype,jtype,:) )
          me%lambda(i,j,:) = lambda
          me%lambda(j,i,:) = lambda
          me%lambda1D(symm1D(i,j),:) = lambda
        end associate
      end do
    end do

    ! Save multithreading-related variables:
    me%nthreads = nthreads
    me%threadCharges = (me%charge%nitems + nthreads - 1)/nthreads
    me%threadTypePairs = (me%nTypePairs + nthreads - 1)/nthreads
    me%threadRigidPairs = (me%nRigidPairs + me%nthreads - 1)/me%nthreads

    ! Compute kspace-model parameters which depend on box geometry:
    call me % update( L )

  end subroutine cKspaceModel_initialize

!---------------------------------------------------------------------------------------------------

  subroutine cKspaceModel_discount_rigid_pairs( me, layer, L, R, Energy, Forces )
    class(cKspaceModel), intent(inout) :: me
    integer,             intent(in)    :: layer
    real(rb),            intent(in)    :: L(3), R(:,:)
    real(rb), optional,  intent(inout) :: Energy, Forces(:,:)

    real(rb) :: E(me%nthreads)

    !$omp parallel num_threads(me%nthreads)
    call discount( omp_get_thread_num() + 1 )
    !$omp end parallel

    if (present(Energy)) Energy = Energy + sum(E)

    contains

      subroutine discount( thread )
        integer, intent(in) :: thread

        integer  :: first, last, k, m
        real(rb) :: invL(3), Rij(3), Fij(3)

        first = (thread - 1)*me%threadTypePairs + 1
        last = min(thread*me%threadTypePairs,me%nTypePairs)

        if (present(Energy)) then
          E(thread) = sum(me%lambda1D(first:last,layer)*me%Erigid(first:last))
        end if

        if (present(Forces)) then
          invL = one/L
          first = (thread - 1)*me%threadRigidPairs + 1
          last = min( thread*me%threadRigidPairs, me%nRigidPairs )
          associate ( lambda => me%lambda1D(:,layer) )
            do k = first, last
              associate ( pair => me%rigidPair(k) )
                Rij = R(:,pair%i) - R(:,pair%j)
                Rij = Rij - L*anint(invL*Rij)
                Fij = lambda(symm1D(pair%itype,pair%jtype))*pair%FbyR*Rij
                do m = 1, 3
                  !$omp atomic update
                  Forces(m,pair%i) = Forces(m,pair%i) + Fij(m)
                  !$omp atomic update
                  Forces(m,pair%j) = Forces(m,pair%j) - Fij(m)
                end do
              end associate
            end do
          end associate
        end if

      end subroutine discount

  end subroutine cKspaceModel_discount_rigid_pairs

!---------------------------------------------------------------------------------------------------

  elemental subroutine cKspaceModel_discount( me, E, W, rsq, QiQj )
    class(cKspaceModel), intent(in)  :: me
    real(rb),            intent(out) :: E, W
    real(rb),            intent(in)  :: rsq, QiQj

    real(rb) :: r, x, expmx2

    r = sqrt(rsq)
    x = me%alpha*r
    expmx2 = exp(-x*x)
    E = -QiQj*uerf( x, expmx2 )/r
    W = E + QiQj*me%beta*expmx2

  end subroutine cKspaceModel_discount

!---------------------------------------------------------------------------------------------------

end module kspaceModelClass
