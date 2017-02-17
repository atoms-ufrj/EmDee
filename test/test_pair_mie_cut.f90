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

#include "emdee.f03"

program test

use EmDee
use mConfig

implicit none

#include "common/declarations.f90"

call command_line_arguments( filename, threads )
call read_data( filename )
call read_configuration( configFile )
call unit_conversions

md = EmDee_system( threads, 1, Rc, Rs, N, c_loc(atomType), c_loc(mass), c_null_ptr )
do i = 1, ntypes
  pair = EmDee_pair_mie_cut( epsilon(i), sigma(i), 12.0_rb, 6.0_rb ) 
  call EmDee_set_pair_model( md, i, i, pair, kCoul )
end do
call EmDee_upload( md, "charges"//c_null_char, c_loc(Q) )
call EmDee_upload( md, "box"//c_null_char, c_loc(L) )
call EmDee_upload( md, "coordinates"//c_null_char, c_loc(R(1,1)) )

call run( 0, Nprop )

contains
  include "common/contained.f90"
end program test

