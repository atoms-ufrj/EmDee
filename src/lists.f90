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
  integer(ib) :: count
  type(c_ptr) :: first
  type(c_ptr) :: last
  type(c_ptr) :: item
end type tList

contains

!---------------------------------------------------------------------------------------------------

  subroutine allocate_list( list, nitems, natoms )
    type(tList), intent(inout) :: list
    integer(ib), intent(in)    :: nitems, natoms

    list%nitems = nitems
    list%count  = 0
    list%item   = malloc_int( nitems )
    list%first  = malloc_int( natoms, value = 1_ib )
    list%last   = malloc_int( natoms, value = 0_ib )

  end subroutine allocate_list

!---------------------------------------------------------------------------------------------------

  subroutine reallocate_list( list, size )
    type(tList), intent(inout) :: list
    integer(ib), intent(in)    :: size

    integer(ib) :: n
    integer(ib), pointer :: old(:), new(:)

    call c_f_pointer( list%item, old, [list%nitems] )
    allocate( new(size) )
    n = min(list%nitems,size)
    new(1:n) = old(1:n)
    deallocate( old )
    list%item = c_loc(new(1))
    list%nitems = size

  end subroutine reallocate_list

!---------------------------------------------------------------------------------------------------

end module lists
