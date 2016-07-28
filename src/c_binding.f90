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

module c_binding

use global

implicit none

contains

!---------------------------------------------------------------------------------------------------

  type(c_ptr) function malloc_int( n, value, array )
    integer(c_int), intent(in)           :: n
    integer(c_int), intent(in), optional :: value, array(:)
    integer(c_int), pointer :: ptr(:)
    allocate( ptr(n) )
    malloc_int = c_loc(ptr(1))
    if (present(array)) then
      ptr = array
    else if (present(value)) then
      ptr = value
    end if
    nullify( ptr )
  end function malloc_int

!---------------------------------------------------------------------------------------------------

  type(c_ptr) function malloc_real( n, value, array )
    integer(c_int), intent(in)           :: n
    real(c_double),    intent(in), optional :: value, array(:)
    real(c_double), pointer :: ptr(:)
    allocate( ptr(n) )
    malloc_real = c_loc(ptr(1))
    if (present(array)) then
      ptr = array
    else if (present(value)) then
      ptr = value
    end if
    nullify( ptr )
  end function malloc_real

!---------------------------------------------------------------------------------------------------

  subroutine realloc_int( ptr, size, new_size )
    type(c_ptr),    intent(inout) :: ptr
    integer(c_int), intent(inout) :: size
    integer(c_int), intent(in)    :: new_size

    integer(c_int) :: n
    integer(c_int), pointer :: old(:), new(:)

    call c_f_pointer( ptr, old, [size] )
    allocate( new(new_size) )
    n = min(size,new_size)
    new(1:n) = old(1:n)
    deallocate( old )
    ptr = c_loc(new(1))
    size = new_size

  end subroutine realloc_int

!---------------------------------------------------------------------------------------------------

  subroutine copy_real( from, to, n, first, last )
    type(c_ptr), intent(in) :: from, to
    integer,     intent(in) :: n, first, last
    real(c_double), pointer :: F(:), T(:)
    call c_f_pointer( from, F, [n] )
    call c_f_pointer( to, T, [n] )
    T(first:last) = F(first:last)
  end subroutine copy_real

!---------------------------------------------------------------------------------------------------

end module c_binding
