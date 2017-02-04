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

module bondModelClass

use global
use modelClass

implicit none

!> An abstract class for bond interaction models:
type, abstract, extends(cModel) :: cBondModel
  contains
    procedure(cBondModel_compute), deferred :: compute
end type cBondModel

!> A class for no-interaction bond model:
type, extends(cBondModel) :: bond_none
  contains
    procedure :: setup => bond_none_setup
    procedure :: compute => bond_none_compute
end type bond_none

abstract interface

  subroutine cBondModel_compute( model, E, W, invR2 )
    import
    class(cBondModel), intent(in)  :: model
    real(rb),          intent(out) :: E, W
    real(rb),          intent(in)  :: invR2
  end subroutine cBondModel_compute

end interface

contains

!===================================================================================================
!                                   B O N D     N O N E
!===================================================================================================

  type(c_ptr) function EmDee_bond_none() bind(C,name="EmDee_bond_none")
    type(bond_none), pointer :: model
    allocate(model)
    call model % setup()
    EmDee_bond_none = model % deliver()
  end function EmDee_bond_none

!---------------------------------------------------------------------------------------------------

  subroutine bond_none_setup( model, params, iparams )
    class(bond_none), intent(inout) :: model
    real(rb), intent(in), optional :: params(:)
    integer,  intent(in), optional :: iparams(:)
    model%name = "none"
  end subroutine bond_none_setup

!---------------------------------------------------------------------------------------------------

  subroutine bond_none_compute( model, E, W, invR2 )
    class(bond_none), intent(in)  :: model
    real(rb),         intent(out) :: E, W
    real(rb),         intent(in)  :: invR2
    E = zero
    W = zero
  end subroutine bond_none_compute

!---------------------------------------------------------------------------------------------------

end module bondModelClass
