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
  character(20) :: name !! TRANSFER TO BASE CLASS
  contains
    procedure(cBondModel_setup), deferred :: setup
    procedure(cBondModel_compute), deferred :: compute
end type cBondModel

!> A class for no-interaction bond model:
type, extends(cBondModel) :: bond_none
  contains
    procedure :: setup => bond_none_setup
    procedure :: compute => bond_none_compute
end type bond_none

!> A container structure for bond models:
type bondModelContainer
  class(cBondModel), allocatable :: model
  contains
    procedure :: bondModelContainer_assign
    generic :: assignment(=) => bondModelContainer_assign
end type bondModelContainer

abstract interface

  subroutine cBondModel_setup( model, params )
    import
    class(cBondModel), intent(inout) :: model
    real(rb),          intent(in)    :: params(:)
  end subroutine cBondModel_setup

  subroutine cBondModel_compute( model, E, F, invR2 )
    import
    class(cBondModel), intent(in)  :: model
    real(rb),          intent(out) :: E, F
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
    call model% setup( [zero] )
    EmDee_bond_none = model % deliver()
  end function EmDee_bond_none

!---------------------------------------------------------------------------------------------------

  subroutine bond_none_setup( model, params )
    class(bond_none), intent(inout) :: model
    real(rb),         intent(in)    :: params(:)
  end subroutine bond_none_setup

!---------------------------------------------------------------------------------------------------

  subroutine bond_none_compute( model, E, F, invR2 )
    class(bond_none), intent(in)  :: model
    real(rb),         intent(out) :: E, F
    real(rb),         intent(in)  :: invR2
    E = zero
    F = zero
  end subroutine bond_none_compute

!===================================================================================================
!                         B O N D     M O D E L    C O N T A I N E R
!===================================================================================================

  subroutine bondModelContainer_assign( new, old )
    class(bondModelContainer), intent(inout) :: new
    type(modelContainer),      intent(in)    :: old

    if (allocated(new%model)) deallocate( new%model )

    if (allocated(old%model)) then
      select type (model => old%model)
        class is (cBondModel)
          allocate( new%model, source = model )
        class default
          stop "ERROR: cannot assign a bond model type from another model type"
      end select
    end if

  end subroutine bondModelContainer_assign

!---------------------------------------------------------------------------------------------------

end module bondModelClass
