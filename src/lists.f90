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

type tList
  integer(ib) :: nitems
  integer(ib) :: nobjects
  integer(ib) :: count
  integer(ib), allocatable :: first(:)
  integer(ib), allocatable :: last(:)
  integer(ib), allocatable :: item(:)
  contains
    procedure :: allocate => tList_allocate
    procedure :: resize => tList_resize
end type tList

contains

!---------------------------------------------------------------------------------------------------

  subroutine tList_allocate( list, nitems, nobjects )
    class(tList), intent(inout) :: list
    integer(ib),  intent(in)    :: nitems, nobjects

    list%nobjects = nobjects
    list%nitems = nitems
    list%count  = 0
    if (allocated(list%first)) deallocate( list%first, list%last, list%item )
    allocate( list%first(nobjects), list%last(nobjects), list%item(nitems) )
    list%first = 1_ib
    list%last = 0_ib

  end subroutine tList_allocate

!---------------------------------------------------------------------------------------------------

  subroutine tList_resize( list, size )
    class(tList), intent(inout) :: list
    integer(ib),  intent(in)    :: size

    integer(ib) :: n
    integer(ib), allocatable :: new(:)

    allocate( new(size) )
    n = min(list%nitems,size)
    new(1:n) = list%item(1:n)
    deallocate( list%item )
    list%item = new
    list%nitems = size

  end subroutine tList_resize

!---------------------------------------------------------------------------------------------------

end module lists
