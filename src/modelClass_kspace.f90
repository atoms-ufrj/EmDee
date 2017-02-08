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
  integer  :: nthreads
  integer  :: threadAtoms

  integer :: natoms, ntypes
  integer,  allocatable :: index(:), type(:)
  real(rb), allocatable :: Q(:)

  integer :: npairs
  integer :: threadPairs
  integer,  allocatable :: pair(:,:)
  real(rb), allocatable :: Edisc(:,:), Wdisc(:,:), FbyR(:)

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

  subroutine cKspaceModel_compute( me, thread, lambda, Potential, Virial, F )
    import
    class(cKspaceModel), intent(inout) :: me
    integer,             intent(in)    :: thread
    real(rb),            intent(in)    :: lambda(:,:)
    real(rb),            intent(inout) :: Potential, Virial
    real(rb), optional,  intent(inout) :: F(:,:)
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
    me%type = types(me%index)
    me%ntypes = maxval(me%type)
    me%Q = charge(me%index)
    me%threadAtoms = (me%natoms + nthreads - 1)/nthreads
    call me % set_parameters( Rc, L, me%alpha )
    call me % update( L )

  end subroutine cKspaceModel_initialize

!---------------------------------------------------------------------------------------------------

  subroutine cKspaceModel_setup_rigid_pairs( me, groupSize, atomIndex, R, charged )
    class(cKspaceModel), intent(inout) :: me
    integer,             intent(in)    :: groupSize(:), atomIndex(:)
    real(rb),            intent(in)    :: R(:,:)
    logical,             intent(in)    :: charged(:)

    integer  :: maxnpairs, npairs, ii, i, jj, j, k, prev, n
    integer  :: ipos, jpos, itype, jtype, imin, imax
    real(rb) :: Ri(3), Qi, rsq, Eij, Wij

    integer,  allocatable :: pos(:), pair(:,:)
    real(rb), allocatable :: FbyR(:)

    allocate( pos(size(charged)) )
    pos(me%index) = [(i,i=1,me%natoms)]

    maxnpairs = sum(groupSize*(groupSize-1)/2)
    allocate( pair(2,maxnpairs), FbyR(maxnpairs) )
    allocate( me%Edisc(me%ntypes,me%ntypes), me%Wdisc(me%ntypes,me%ntypes) )

    npairs = 0
    prev = 0
    me%Edisc = zero
    me%Wdisc = zero
    do k = 1, size(groupSize)
      n = groupSize(k)
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
              imin = min(itype,jtype)
              imax = max(itype,jtype)

              rsq = sum((Ri - R(:,j))**2)
              call me % discount( Eij, Wij, rsq, Qi*me%Q(jpos) )
              me%Edisc(imin,imax) = me%Edisc(imin,imax) + Eij
              me%Wdisc(imin,imax) = me%Wdisc(imin,imax) + Wij

              npairs = npairs + 1
              pair(:,npairs) = [ipos, jpos]
              FbyR(npairs) = Wij/rsq
            end if
          end do
        end if
      end do
      prev = prev + n
    end do

    me%npairs = npairs
    me%threadPairs = (npairs + me%nthreads - 1)/me%nthreads
    me%pair = pair(:,1:npairs)
    me%FbyR = FbyR(1:npairs)

  end subroutine cKspaceModel_setup_rigid_pairs

!---------------------------------------------------------------------------------------------------

  subroutine cKspaceModel_discount_rigid_pairs( me, thread, lambda, R, Potential, Virial, F )
    class(cKspaceModel), intent(inout) :: me
    integer,             intent(in)    :: thread
    real(rb),            intent(in)    :: lambda(:,:), R(:,:)
    real(rb),            intent(inout) :: Potential, Virial
    real(rb), optional,  intent(inout) :: F(:,:)

    integer  :: i, j, k, ii, jj, itype, jtype
    real(rb) :: Fij(3)

    if (present(F)) then
      do k = (thread - 1)*me%threadPairs + 1, min(thread*me%threadPairs,me%npairs)
        ii = me%pair(1,k)
        jj = me%pair(2,k)
        itype = me%type(ii)
        jtype = me%type(jj)
        i = me%index(ii)
        j = me%index(jj)
        Fij = lambda(itype,jtype)*me%FbyR(k)*(R(:,i) - R(:,j))
        F(:,i) = F(:,i) + Fij
        F(:,j) = F(:,j) - Fij
      end do
    end if

    if (thread == 1) then
      do itype = 1, me%ntypes
        do jtype = itype, me%ntypes
          Potential = Potential + lambda(itype,jtype)*me%Edisc(itype,jtype)
          Virial = Virial + lambda(itype,jtype)*me%Wdisc(itype,jtype)
        end do
      end do
    end if

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
    E = -QiQj*uerf( x, expmx2 )
    W = QiQj*prefactor*x*expmx2

  end subroutine cKspaceModel_discount

!---------------------------------------------------------------------------------------------------

end module kspaceModelClass
