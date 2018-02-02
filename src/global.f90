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

module global

use, intrinsic :: iso_c_binding

implicit none

integer, parameter :: ib = c_int
integer, parameter :: rb = c_double
integer, parameter :: lb = c_bool
integer, parameter :: sl = 256

real(rb), parameter :: zero   = 0.0_rb,                 &
                       one    = 1.0_rb,                 &
                       two    = 2.0_rb,                 &
                       three  = 3.0_rb,                 &
                       half   = 0.5_rb,                 &
                       third  = one/three,              &
                       fourth = 0.25_rb,                &
                       sixth  = 1.0_rb/6.0_rb,          &
                       pi     = 3.14159265358979324_rb, &
                       twoPi  = two*pi,                 &
                       piBy2  = 0.5_rb*pi,              &
                       invSqrt3 = one/sqrt(3.0_rb)

contains

!===================================================================================================

  subroutine error( task, msg )
    use, intrinsic :: iso_fortran_env
    character(*), intent(in) :: task, msg
    write(ERROR_UNIT,'("Error in ",A,": ",A,".")') trim(task), trim(msg)
    call exit(1)
  end subroutine error

!===================================================================================================

  subroutine warning( msg )
    use, intrinsic :: iso_fortran_env
    character(*), intent(in) :: msg
    write(ERROR_UNIT,'("WARNING: ",A,".")') trim(msg)
  end subroutine warning

!===================================================================================================

  function ranged( i, imax )
    integer, intent(in) :: i(:), imax
    logical             :: ranged
    integer :: j
    j = 0
    ranged = .true.
    do while (ranged.and.(j < size(i)))
      j = j + 1
      ranged = (i(j) > 0).and.(i(j) <= imax)
    end do
  end function ranged

!===================================================================================================

  function sorted( x, indices ) result( y )
    integer, intent(in)           :: x(:)
    logical, intent(in), optional :: indices
    integer                       :: y(size(x))
    integer :: a(size(x)), ind(size(x)), N, i, j
    logical :: swapped
    N = size(x)
    a = x
    ind = [(i,i=1,N)]
    do j = N-1, 1, -1
      swapped = .false.
      do i = 1, j
        if (a(i) > a(i+1)) then ! increasing order
!        if (a(i) < a(i+1)) then ! decreasing order
          call swap( a(i), a(i+1) )
          call swap( ind(i), ind(i+1) )
          swapped = .true.
        end if
      end do
      if (.not. swapped) exit
    end do
    if (present(indices)) then
      y = merge(ind,a,indices)
    else
      y = a
    end if
    contains
      subroutine swap( a, b )
        integer, intent(inout) :: a, b
        integer :: temp
        temp = a
        a = b
        b = temp
      end subroutine swap
  end function sorted

!===================================================================================================

  character(sl) function string( carray )
    character(c_char), intent(in) :: carray(*)
    integer :: i
    string = ""
    do i = 1, sl
      if (carray(i) == c_null_char) return
      string(i:i) = carray(i)
    end do
  end function string

!===================================================================================================

end module global
