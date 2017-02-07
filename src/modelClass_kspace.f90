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

implicit none

!> An abstract class for kspace interaction models:
type, abstract, extends(cModel) :: cKspaceModel

  real(rb) :: alpha
  integer  :: nthreads
  integer  :: threadAtoms

  integer :: natoms, ntypes
  integer,  allocatable :: index(:), type(:)
  real(rb), allocatable :: Q(:)

  contains
    procedure :: initialize => cKspaceModel_initialize
    procedure(cKspaceModel_set_parameters), deferred :: set_parameters
    procedure(cKspaceModel_update),  deferred :: update
    procedure(cKspaceModel_prepare), deferred :: prepare
    procedure(cKspaceModel_compute), deferred :: compute
end type cKspaceModel

abstract interface

  subroutine cKspaceModel_set_parameters( me, Rc, L )
    import
    class(cKspaceModel), intent(inout) :: me
    real(rb),            intent(in)    :: Rc, L(3)
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

  subroutine cKspaceModel_compute( me, thread, lambda, F, Potential, Virial )
    import
    class(cKspaceModel), intent(inout) :: me
    integer,             intent(in)    :: thread
    real(rb),            intent(in)    :: lambda(:,:)
    real(rb),            intent(inout) :: F(:,:), Potential, Virial
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
    call me % set_parameters( Rc, L )
    call me % update( L )

  end subroutine cKspaceModel_initialize

!---------------------------------------------------------------------------------------------------

end module kspaceModelClass
