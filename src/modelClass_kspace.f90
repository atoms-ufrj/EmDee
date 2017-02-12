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
use modelClass
use omp_lib
use math

implicit none

real(rb), parameter, private :: prefactor = two/sqrt(Pi)

!> An abstract class for kspace interaction models:
type, abstract, extends(cModel) :: cKspaceModel

  real(rb) :: alpha
  real(rb) :: beta

  integer :: natoms                       ! Number of charged atoms
  integer,  allocatable :: index(:)       ! System-wise index of each charged atom
  integer,  allocatable :: type(:)        ! Type of each charged atom
  real(rb), allocatable :: Q(:)           ! Charge of each charged atom

  integer :: ntypes                       ! Number of types of charged atoms
  integer :: npairtypes                   ! Number of types of atoms pairs
  integer :: npairs                       ! Number of intrabody pairs (fixed distance)
  integer,  allocatable :: pair(:,:)      ! Local indices of atoms forming each intrabody pair
  real(rb), allocatable :: FbyR(:)        ! Unscaled force/distance ratio for each intrabody pair
  real(rb), allocatable :: Erigid(:)      ! Unscaled energy for intrabody pairs of each pair type
  real(rb), allocatable :: Wrigid(:)      ! Unscaled virial for intrabody pairs of each pair type

  integer :: nthreads
  integer :: threadAtoms
  integer :: threadPairs
  integer :: threadPairTypes

  contains
    procedure :: initialize => cKspaceModel_initialize
    procedure :: setup_rigid_pairs => cKspaceModel_setup_rigid_pairs
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

  subroutine cKspaceModel_prepare( me, thread, Rs )
    import
    class(cKspaceModel), intent(inout) :: me
    integer,             intent(in)    :: thread
    real(rb),            intent(in)    :: Rs(:,:)
  end subroutine cKspaceModel_prepare

  subroutine cKspaceModel_compute( me, thread, compute, lambda, E, W, F )
    import
    class(cKspaceModel), intent(inout) :: me
    integer(ib),         intent(in)    :: thread
    logical(lb),         intent(in)    :: compute
    real(rb),            intent(in)    :: lambda(:,:)
    real(rb),            intent(out)   :: E
    real(rb),            intent(inout) :: W, F(:,:)
  end subroutine cKspaceModel_compute

end interface

contains

!---------------------------------------------------------------------------------------------------

  subroutine cKspaceModel_initialize( me, nthreads, Rc, L, types, charged, charge )
    class(cKspaceModel), intent(inout) :: me
    integer,             intent(in)    :: nthreads
    real(rb),            intent(in)    :: Rc, L(3)
    integer,             intent(in)    :: types(:)
    logical,             intent(in)    :: charged(size(types))
    real(rb),            intent(in)    :: charge(size(types))

    character(*), parameter :: task = "kspace model initialization"

    integer :: i

    me%nthreads = nthreads
    me%index = pack( [(i,i=1,size(types))], charged )
    me%natoms = size(me%index)
    if (me%natoms == 0) call error( task, "system without charged atoms" )
    me%threadAtoms = (me%natoms + nthreads - 1)/nthreads
    me%type = types(me%index)
    me%ntypes = maxval(me%type)
    me%Q = charge(me%index)

    call me % set_parameters( Rc, L, me%alpha )
    me%beta = two*me%alpha/sqrt(Pi)

    call me % update( L )

    me%npairtypes = me%ntypes*(me%ntypes - 1)/2 + me%ntypes
    me%threadPairTypes = (me%npairtypes + nthreads - 1)/nthreads
    allocate( me%Erigid(me%npairtypes), me%Wrigid(me%npairtypes) )

  end subroutine cKspaceModel_initialize

!---------------------------------------------------------------------------------------------------

  subroutine cKspaceModel_setup_rigid_pairs( me, groupSize, atomIndex, R, charged )
    class(cKspaceModel), intent(inout) :: me
    integer,             intent(in)    :: groupSize(:), atomIndex(:)
    real(rb),            intent(in)    :: R(:,:)
    logical,             intent(in)    :: charged(:)

    integer  :: maxnpairs, npairs, ii, i, jj, j, kk, k, prev, n
    integer  :: ipos, jpos, itype, jtype, pos(me%natoms)
    real(rb) :: Ri(3), Qi, rsq, Eij, Wij, Eself(me%ntypes)

    integer,  allocatable :: pair(:,:)
    real(rb), allocatable :: FbyR(:)

    pos = [(i,i=1,me%natoms)]
    pos(me%index) = pos

    maxnpairs = sum(groupSize*(groupSize-1)/2)
    allocate( pair(2,maxnpairs), FbyR(maxnpairs) )

    npairs = 0
    prev = 0
    me%Erigid = zero
    me%Wrigid = zero
    do kk = 1, size(groupSize)
      n = groupSize(kk)
      do ii = 1, n-1
        i = atomIndex(prev+ii)
        if (charged(i)) then
          Ri = R(:,i)
          ipos = pos(i)
          Qi = me%Q(ipos)
          itype = me%type(ipos)
          do jj = ii+1, n
            j = atomIndex(prev+jj)
            if (charged(j)) then
              jpos = pos(j)
              jtype = me%type(jpos)

              rsq = sum((Ri - R(:,j))**2)
              call me % discount( Eij, Wij, rsq, Qi*me%Q(jpos) )

              k = symm1D( itype, jtype )
              me%Erigid(k) = me%Erigid(k) + Eij
              me%Wrigid(k) = me%Wrigid(k) + Wij

              npairs = npairs + 1
              pair(:,npairs) = [ipos, jpos]
              FbyR(npairs) = Wij/rsq
            end if
          end do
        end if
      end do
      prev = prev + n
    end do

    Eself = zero
    do i = 1, me%natoms
      itype = me%type(i)
      Eself(itype) = Eself(itype) + me%Q(i)**2
    end do
    Eself = -me%alpha/sqrt(Pi)*Eself
    do itype = 1, me%ntypes
      k = symm1D( itype, itype )
      me%Erigid(k) = me%Erigid(k) + Eself(itype)
    end do

    me%npairs = npairs
    me%threadPairs = (npairs + me%nthreads - 1)/me%nthreads
    me%pair = pair(:,1:npairs)
    me%FbyR = FbyR(1:npairs)

  end subroutine cKspaceModel_setup_rigid_pairs

!---------------------------------------------------------------------------------------------------

  subroutine cKspaceModel_discount_rigid_pairs( me, thread, lambda, R, F )
    class(cKspaceModel), intent(inout) :: me
    integer,             intent(in)    :: thread
    real(rb),            intent(in)    :: lambda(:), R(:,:)
    real(rb),            intent(inout) :: F(:,:)

    integer  :: i, j, k, m, ii, jj
    real(rb) :: Fij(3)

    do k = (thread - 1)*me%threadPairs + 1, min(thread*me%threadPairs,me%npairs)
      ii = me%pair(1,k)
      jj = me%pair(2,k)
      i = me%index(ii)
      j = me%index(jj)
      Fij = lambda(symm1D(me%type(ii),me%type(ii)))*me%FbyR(k)*(R(:,i) - R(:,j))
      do m = 1, 3
        !$omp atomic update
        F(m,i) = F(m,i) + Fij(m)
        !$omp atomic update
        F(m,j) = F(m,j) - Fij(m)
      end do
    end do

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
