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

program test

use iso_c_binding

implicit none

#include "emdee.f03"

integer, parameter :: ib = 4, rb = 8

integer(ib) :: N, Nsteps, Nprop
real(rb)    :: rho, Rc, Rs, Rc2, Temp, Dt, Dt_2
real(rb), target  :: L
real(rb), pointer :: R(:,:), V(:,:), Q(:)
integer,  pointer :: types(:)

type(tEmDee), target :: md
type(c_ptr),  target :: pair

integer :: threads
character(256) :: filename

call command_line_arguments( filename, threads )
call read_data( filename )
call create_configuration
call set_charges( types, Q )

md = EmDee_system( threads, 1, Rc, Rs, N, c_loc(types), c_null_ptr )
pair = EmDee_pair_softcore_cut( 1.0_rb, 1.0_rb, 1.0_rb )

call EmDee_set_pair_model( md, 1, 1, pair )
call EmDee_set_pair_model( md, 2, 2, pair )

call EmDee_upload( md, "charges"//c_null_char, c_loc(Q) )
call EmDee_upload( md, "box"//c_null_char, c_loc(L) )
call EmDee_upload( md, "coordinates"//c_null_char, c_loc(R(1,1)) )

call run( 0, Nprop )

contains
  include "common_routines.inc"
end program test

