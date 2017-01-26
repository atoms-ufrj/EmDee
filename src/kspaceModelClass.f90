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

implicit none

!> An abstract class for kspace interaction models:
type, abstract, extends(cModel) :: cKspaceModel

  real(rb) :: alpha
  integer  :: nthreads
  logical  :: verbose

  real(rb) :: E_self

  integer :: natoms
  integer,  allocatable :: index(:)
  real(rb), allocatable :: Q(:)

  integer :: atoms_per_thread

  contains
    procedure :: initialize => cKspaceModel_initialize
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

    model%index = pack( [(i,i=1,size(Q))], Q /= zero )
    model%natoms = size(model%index)
    if (model%natoms == 0) call error( "kspace model initialization", "no charged atoms" )
    model%Q = Q(model%index)

    model%atoms_per_thread = (model%natoms + nthreads - 1)/nthreads

    call model % set_parameters( Rc, L, Q )

    call model % update( L )

    model%E_self = model%alpha*sum(model%Q**2)/sqrt(Pi)

  end subroutine cKspaceModel_initialize

!---------------------------------------------------------------------------------------------------

end module kspaceModelClass
