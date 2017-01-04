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
  integer :: nthreads
  contains
    procedure(cKspaceModel_initialize), deferred :: initialize
    procedure(cKspaceModel_compute), deferred :: compute
end type cKspaceModel

abstract interface

  subroutine cKspaceModel_initialize( model, nthreads, Rc, invL, q )
    import
    class(cKspaceModel), intent(inout) :: model
    integer,             intent(in)    :: nthreads
    real(rb),            intent(in)    :: Rc, invL(3), q(:)
  end subroutine cKspaceModel_initialize

  subroutine cKspaceModel_compute( me, nthreads, q, Rs, V, E )
    import
    class(cKspaceModel), intent(in)  :: me
    integer,             intent(in)  :: nthreads
    real(rb),            intent(in)  :: q(:), Rs(3,size(q)), V
    real(rb),            intent(out) :: E
  end subroutine cKspaceModel_compute

end interface

!> A class for no-interaction kspace model:
!type, extends(cKspaceModel) :: kspace_none
!  contains
!    procedure :: setup => kspace_none_setup
!    procedure :: compute => kspace_none_compute
!end type kspace_none

contains

!===================================================================================================
!                                   K S P A C E     N O N E
!===================================================================================================

!  type(c_ptr) function EmDee_kspace_none() bind(C,name="EmDee_kspace_none")
!    type(kspace_none), pointer :: model
!    allocate(model)
!    call model % setup()
!    EmDee_kspace_none = model % deliver()
!  end function EmDee_kspace_none

!!---------------------------------------------------------------------------------------------------

!  subroutine kspace_none_setup( model, params, iparams )
!    class(kspace_none), intent(inout) :: model
!    real(rb), intent(in), optional :: params(:)
!    integer,  intent(in), optional :: iparams(:)
!    model%name = "none"
!  end subroutine kspace_none_setup

!!---------------------------------------------------------------------------------------------------

!  subroutine kspace_none_compute( me, nthreads, q, Rs, V, E )
!    class(kspace_none), intent(in)  :: me
!    integer,            intent(in)  :: nthreads
!    real(rb),           intent(in)  :: q(:), Rs(3,size(q)), V
!    real(rb),           intent(out) :: E
!    E = zero
!  end subroutine kspace_none_compute

!---------------------------------------------------------------------------------------------------

end module kspaceModelClass
