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

module lists

use global

implicit none

type tList
  integer :: nitems
  integer :: nobjects
  integer :: count
  integer,  pointer     :: first(:) => NULL()
  integer,  pointer     :: middle(:) => NULL()
  integer,  pointer     :: last(:) => NULL()
  integer,  allocatable :: item(:)
  logical,  allocatable :: check(:)
  real(rb), allocatable :: value(:)
  contains
    procedure :: allocate => tList_allocate
    procedure :: resize => tList_resize
    procedure :: object => tList_object
end type tList

contains

!---------------------------------------------------------------------------------------------------

  elemental subroutine tList_allocate( list, nitems, nobjects, middle, check, value )
    class(tList),      intent(inout) :: list
    integer,           intent(in)    :: nitems, nobjects
    logical, optional, intent(in)    :: middle, check, value

    list%nobjects = nobjects
    list%nitems = nitems
    list%count  = 0
    if (associated(list%first)) deallocate( list%first, list%last, list%item )
    allocate( list%first(nobjects), list%last(nobjects), list%item(nitems) )
    list%first = 1
    list%last = 0
    if (present(middle)) then
      if (middle) allocate( list%middle(nobjects) )
    end if
    if (present(check)) then
      if (check) allocate( list%check(nitems) )
    end if
    if (present(value)) then
      if (value) allocate( list%value(nitems) )
    end if
  end subroutine tList_allocate

!---------------------------------------------------------------------------------------------------

  subroutine tList_resize( list, size )
    class(tList), intent(inout) :: list
    integer,      intent(in)    :: size

    integer :: n
    integer,  allocatable :: item(:)
    logical,  allocatable :: check(:)
    real(rb), allocatable :: value(:)

    allocate( item(size) )
    n = min(list%nitems,size)
    item(1:n) = list%item(1:n)
    deallocate( list%item )
    call move_alloc( item, list%item )
    list%nitems = size
    if (allocated(list%check)) then
      allocate( check(size) )
      check(1:n) = list%check(1:n)
      deallocate( list%check )
      call move_alloc( check, list%check )
    end if
    if (allocated(list%value)) then
      allocate( value(size) )
      value(1:n) = list%value(1:n)
      deallocate( list%value )
      call move_alloc( value, list%value )
    end if

  end subroutine tList_resize

!---------------------------------------------------------------------------------------------------

  pure function tList_object( list, item ) result( object )
    class(tList), intent(in) :: list
    integer,      intent(in) :: item
    integer                  :: object
    object = 1
    do while (list%last(object) < item)
      object = object + 1
    end do
  end function tList_object

!---------------------------------------------------------------------------------------------------

end module lists
