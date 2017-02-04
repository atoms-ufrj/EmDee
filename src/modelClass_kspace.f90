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
use lists
use omp_lib

implicit none

!> An abstract class for kspace interaction models:
type, abstract, extends(cModel) :: cKspaceModel

  real(rb) :: alpha
  integer  :: nthreads
  logical  :: verbose

  real(rb) :: beta
  real(rb) :: E_self

  integer :: natoms
  integer,  allocatable :: index(:)
  logical,  allocatable :: charged(:)
  real(rb), allocatable :: Q(:)

  integer :: atoms_per_thread

  contains
    procedure :: initialize => cKspaceModel_initialize
    procedure :: discount => cKspaceModel_discount
    procedure(cKspaceModel_set_parameters), deferred :: set_parameters
    procedure(cKspaceModel_update), deferred :: update
    procedure(cKspaceModel_compute), deferred :: compute
end type cKspaceModel

abstract interface

  subroutine cKspaceModel_set_parameters( model, Rc, L, Q )
    import
    class(cKspaceModel), intent(inout) :: model
    real(rb),            intent(in)    :: Rc, L(3), Q(:)
  end subroutine cKspaceModel_set_parameters

  ! This procedure must be invoked whenever the box geometry or atoms charges change: 
  subroutine cKspaceModel_update( model, L )
    import
    class(cKspaceModel), intent(inout) :: model
    real(rb),            intent(in)    :: L(3)
  end subroutine cKspaceModel_update

  ! This procedure computes the reciprocal part of Ewald-type electrostatics:
  subroutine cKspaceModel_compute( model, R, F, Potential, Virial )
    import
    class(cKspaceModel), intent(in)    :: model
    real(rb),            intent(in)    :: R(:,:)
    real(rb),            intent(inout) :: F(:,:), Potential, Virial
  end subroutine cKspaceModel_compute

end interface

contains

!---------------------------------------------------------------------------------------------------

  subroutine cKspaceModel_initialize( model, nthreads, Rc, L, Q, verbose )
    class(cKspaceModel), intent(inout) :: model
    integer,             intent(in)    :: nthreads
    real(rb),            intent(in)    :: Rc, L(3), Q(:)
    logical,             intent(in)    :: verbose

    integer :: i

    model%nthreads = nthreads
    model%verbose = verbose

    model%charged = Q /= zero
    model%index = pack( [(i,i=1,size(Q))], model%charged )
    model%natoms = size(model%index)
    if (model%natoms == 0) call error( "kspace model initialization", "no charged atoms" )
    model%Q = Q

    model%atoms_per_thread = (model%natoms + nthreads - 1)/nthreads

    call model % set_parameters( Rc, L, Q )
    call model % update( L )

    model%beta = two*model%alpha/sqrt(Pi)
    model%E_self = model%alpha*sum(model%Q**2)/sqrt(Pi)

  end subroutine cKspaceModel_initialize

!---------------------------------------------------------------------------------------------------

  subroutine cKspaceModel_discount( model, thread, excluded, R, F, Potential, Virial )
    class(cKspaceModel), intent(inout) :: model
    integer,             intent(in)    :: thread
    type(tList),         intent(in)    :: excluded
    real(rb),            intent(in)    :: R(:,:)
    real(rb),            intent(inout) :: F(:,:), Potential, Virial

    integer  :: first, last, ii, jj, i, j
    real(rb) :: Qi, Ri(3), Rij(3), R2, Rd, QiQj, alphaRd, Eij, Wij, Fij(3)

    first = (thread - 1)*model%atoms_per_thread + 1
    last = min(thread*model%atoms_per_thread, model%natoms)

    do ii = first, last
      i = model%index(i)
      Qi = model%Q(ii)
      Ri = R(:,i)
      do jj = excluded%first(i), excluded%last(i)
        j = excluded%item(jj)
        if (model%charged(j)) then
          Rij = Ri - R(:,j)
          R2 = sum(Rij*Rij)
          Rd = sqrt(R2)
          QiQj = Qi*model%Q(j)
          alphaRd = model%alpha*Rd
          Eij = QiQj*erf(alphaRd)/Rd
          Wij = Eij - QiQj*model%beta*exp(-alphaRd*alphaRd)
          Fij = Wij*Rij/R2
          F(:,i) = F(:,i) - Fij
          F(:,j) = F(:,j) + Fij
          Potential = Potential - Eij
          Virial = Virial - Wij
        end if
      end do
    end do

 end subroutine cKspaceModel_discount

!---------------------------------------------------------------------------------------------------

end module kspaceModelClass
