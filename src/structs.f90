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

module structs

use models

implicit none

integer(ib), parameter, private :: extra = 500

type tStruct
  integer(ib) :: i, j, k, l
  class(cModel), allocatable :: model
end type tStruct

type structList
  integer :: number = 0
  integer :: max = 0
  logical :: exist = .false.
  type(tStruct), allocatable :: item(:)
  contains
    procedure :: add => structList_add
end type structList

contains

!---------------------------------------------------------------------------------------------------

  subroutine structList_add( struct, i, j, k, l, model )
    class(structList),      intent(inout) :: struct
    integer(ib),            intent(in)    :: i, j, k, l
    class(cModel), pointer, intent(in)    :: model

    type(tStruct), allocatable :: new(:)

    if (.not.struct%exist) then
      struct%max = extra
      allocate( struct%item(struct%max) )
      struct%exist = .true.
    else if (struct%number == struct%max) then
      struct%max = struct%max + extra
      allocate( new(struct%max) )
      new(1:struct%number) = struct%item
      deallocate( struct%item )
      call move_alloc( new, struct%item )
    end if
    struct%number = struct%number + 1
    struct%item(struct%number) = tStruct( i, j, k, l, model )

  end subroutine structList_add

!---------------------------------------------------------------------------------------------------

  subroutine add_bonded_struc( struct, i, j, k, l, model )
    type(c_ptr),            intent(inout) :: struct
    integer(ib),            intent(in)    :: i, j, k, l
    class(cModel), pointer, intent(in)    :: model

    type(structList), pointer :: ptr
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
    ptr%item(ptr%number)%i = i
    ptr%item(ptr%number)%j = j
    ptr%item(ptr%number)%k = k
    ptr%item(ptr%number)%l = l
    allocate( ptr%item(ptr%number)%model, source = model )
  end subroutine add_bonded_struc

!---------------------------------------------------------------------------------------------------

end module structs
