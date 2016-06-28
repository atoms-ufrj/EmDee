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

module bonded_structs

use c_binding
use models

implicit none

integer(ib), parameter, private :: extra = 500

type tStruct
  integer(ib) :: i, j, k, l
  type(EmDee_Model), pointer :: model => null()
end type tStruct

type tStructData
  integer :: number = 0
  integer :: max = 0
  type(tStruct), allocatable :: item(:)
end type tStructData

contains

!---------------------------------------------------------------------------------------------------

  subroutine add_bonded_struc( struct, i, j, k, l, model )
    type(c_ptr), intent(inout) :: struct
    integer(ib), intent(in)    :: i, j, k, l
    type(c_ptr), intent(in)    :: model

    type(tStructData), pointer :: ptr
    type(tStruct), allocatable :: new(:)

    if (c_associated(struct)) then
      call c_f_pointer( struct, ptr )
    else
      allocate( ptr )
      allocate( ptr%item(0) )
      struct = c_loc(ptr)
    end if

    if (ptr%number == ptr%max) then
      allocate( new(ptr%max+extra) )
      new(1:ptr%number) = ptr%item
      deallocate( ptr%item )
      call move_alloc( new, ptr%item )
    end if
    ptr%number = ptr%number + 1
    ptr%item(ptr%number) = tStruct( i, j, k, l )
    call c_f_pointer( model, ptr%item(ptr%number)%model )
  end subroutine add_bonded_struc

!---------------------------------------------------------------------------------------------------

end module bonded_structs
