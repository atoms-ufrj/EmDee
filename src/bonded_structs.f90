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

type, bind(C) :: tStructData
  integer     :: number
  integer     :: max
  type(c_ptr) :: list
  real(rb)    :: energy
  real(rb)    :: virial
end type tStructData

contains

!---------------------------------------------------------------------------------------------------

  subroutine add_bonded_struc( struct, i, j, k, l, model )
    type(tStructData), intent(inout) :: struct
    integer(ib),       intent(in)    :: i, j, k, l
    type(c_ptr),       intent(in)    :: model

    type(tStruct), pointer :: old(:), new(:)

    if (struct%number + 1 > struct%max) then
      call c_f_pointer( struct%list, old, [struct%max] )
      allocate( new(struct%max+extra) )
      new(1:struct%number) = old(1:struct%number)
      deallocate( old )
      struct%list = c_loc(new(1))
    else
      call c_f_pointer( struct%list, new, [struct%number+1] )
    end if
    struct%number = struct%number + 1
    new(struct%number) = tStruct( i, j, k, l )
    call c_f_pointer( model, new(struct%number)%model )

    nullify( new )
  end subroutine add_bonded_struc

!---------------------------------------------------------------------------------------------------

end module bonded_structs
