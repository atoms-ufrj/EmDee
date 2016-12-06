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

module angleModelClass

use global
use modelClass

implicit none

!> An abstract class for angle interaction models:
type, abstract, extends(cModel) :: cAngleModel
  contains
    procedure(cAngleModel_compute), deferred :: compute
end type cAngleModel

!> A class for no-interaction angle model:
type, extends(cAngleModel) :: angle_none
  contains
    procedure :: setup => angle_none_setup
    procedure :: compute => angle_none_compute
end type angle_none

abstract interface

  subroutine cAngleModel_compute( model, Ea, Fa, theta )
    import
    class(cAngleModel), intent(in) :: model
    real(rb),          intent(out) :: Ea, Fa
    real(rb),          intent(in)  :: theta
  end subroutine cAngleModel_compute

end interface

contains

!===================================================================================================
!                                   A N G L E     N O N E
!===================================================================================================

  type(c_ptr) function EmDee_angle_none() bind(C,name="EmDee_angle_none")
    type(angle_none), pointer :: model
    allocate(model)
    call model % setup()
    EmDee_angle_none = model % deliver()
  end function EmDee_angle_none

!---------------------------------------------------------------------------------------------------

  subroutine angle_none_setup( model, params, iparams )
    class(angle_none), intent(inout) :: model
    real(rb), intent(in), optional :: params(:)
    integer,  intent(in), optional :: iparams(:)
    model%name = "none"
  end subroutine angle_none_setup

!---------------------------------------------------------------------------------------------------

  subroutine angle_none_compute( model, Ea, Fa, theta )
    class(angle_none), intent(in)  :: model
    real(rb),          intent(out) :: Ea, Fa
    real(rb),          intent(in)  :: theta
    Ea = zero
    Fa = zero
  end subroutine angle_none_compute

!---------------------------------------------------------------------------------------------------

end module angleModelClass
