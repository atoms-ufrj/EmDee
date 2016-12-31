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

module modelClass

use global
use, intrinsic :: iso_c_binding

type, abstract :: cModel
  character(20) :: name = "none"
  contains
    procedure :: deliver => cModel_deliver
    procedure(cModel_setup), deferred :: setup
end type cModel

abstract interface

  subroutine cModel_setup( model, params, iparams )
    import
    class(cModel), intent(inout) :: model
    real(rb), intent(in), optional :: params(:)
    integer,  intent(in), optional :: iparams(:)
  end subroutine cModel_setup

end interface

type modelContainer
  class(cModel), allocatable :: model
  contains
    procedure :: modelContainer_assign
    generic :: assignment(=) => modelContainer_assign
end type modelContainer

contains

!===================================================================================================
!                                   M O D E L    C L A S S
!===================================================================================================

  type(c_ptr) function cModel_deliver( this )
    class(cModel), intent(in) :: this

    type(modelContainer), pointer :: container

    allocate( container )
    allocate( container%model, source = this )
    cModel_deliver = c_loc(container)

  end function cModel_deliver

!===================================================================================================
!                                M O D E L    C O N T A I N E R
!===================================================================================================

  subroutine modelContainer_assign( new, old )
    class(modelContainer), intent(inout) :: new
    type(modelContainer),  intent(in)    :: old

    if (allocated(new%model)) deallocate( new%model )

    if (allocated(old%model)) then
      allocate( new%model, source = old%model )
    end if

  end subroutine modelContainer_assign

!---------------------------------------------------------------------------------------------------

end module modelClass
