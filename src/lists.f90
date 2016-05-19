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
!            Federal University of Rio de Janeiro, Brazilmodule lists

module lists

use c_binding

implicit none

type, bind(C) :: tList
  integer(ib) :: nitems
  integer(ib) :: nobjects
  integer(ib) :: count
  type(c_ptr) :: first = c_null_ptr
  type(c_ptr) :: last = c_null_ptr
  type(c_ptr) :: item = c_null_ptr
end type tList

type tListPtr
  integer(ib), pointer :: first(:)
  integer(ib), pointer :: last(:)
  integer(ib), pointer :: item(:)
end type tListPtr

contains

!---------------------------------------------------------------------------------------------------

  subroutine allocate_list( list, nitems, nobjects )
    type(tList), intent(inout) :: list
    integer(ib), intent(in)    :: nitems, nobjects

    type(tListPtr) :: ptr

    list%nobjects = nobjects
    list%nitems = nitems
    list%count  = 0
    if (c_associated(list%first)) then
      call c_f_list( list, ptr )
      deallocate( ptr%first, ptr%last, ptr%item )
    end if
    list%first  = malloc_int( nobjects, value = 1_ib )
    list%last   = malloc_int( nobjects, value = 0_ib )
    list%item   = malloc_int( nitems )

  end subroutine allocate_list

!---------------------------------------------------------------------------------------------------

  subroutine resize_list( list, ptr, size )
    type(tList),    intent(inout) :: list
    type(tListPtr), intent(inout) :: ptr
    integer(ib),    intent(in)    :: size

    integer(ib) :: n
    integer(ib), pointer :: old(:), new(:)

    call c_f_pointer( list%item, old, [list%nitems] )
    allocate( new(size) )
    n = min(list%nitems,size)
    new(1:n) = old(1:n)
    deallocate( old )
    list%item = c_loc(new(1))
    list%nitems = size
    call c_f_pointer( list%item, ptr%item, [list%nitems] )

  end subroutine resize_list

!---------------------------------------------------------------------------------------------------

  subroutine c_f_list( list, ptr )
    type(tList),    intent(in)  :: list
    type(tListPtr), intent(out) :: ptr
    call c_f_pointer( list%first, ptr%first, [list%nobjects] )
    call c_f_pointer( list%last,  ptr%last,  [list%nobjects] )
    call c_f_pointer( list%item,  ptr%item,  [list%count] )
  end subroutine c_f_list

!---------------------------------------------------------------------------------------------------

end module lists
