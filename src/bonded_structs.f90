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
  type(tModel), pointer :: model => null()
end type tStruct

contains

!---------------------------------------------------------------------------------------------------

  subroutine add_bonded_struc( ptr, size, sizemax, i, j, k, l, model )
    type(c_ptr), intent(inout) :: ptr
    integer(ib), intent(inout) :: size, sizemax
    integer(ib), intent(in)    :: i, j, k, l
    type(c_ptr), intent(in)    :: model

    type(tStruct), pointer :: old(:), new(:)

    if (size + 1 > sizemax) then
      call c_f_pointer( ptr, old, [sizemax] )
      allocate( new(sizemax+extra) )
      new(1:size) = old(1:size)
      deallocate( old )
      ptr = c_loc(new(1))
    else
      call c_f_pointer( ptr, new, [size+1] )
    end if
    size = size + 1
    new(size) = tStruct( i, j, k, l )
    call c_f_pointer( model, new(size)%model )

    nullify( new )
  end subroutine add_bonded_struc

!---------------------------------------------------------------------------------------------------

end module bonded_structs
