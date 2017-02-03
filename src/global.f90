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
                       half   = 0.5_rb,                 &
                       third  = 0.33333333333333333_rb, &
                       fourth = 0.25_rb,                &
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

end module global
