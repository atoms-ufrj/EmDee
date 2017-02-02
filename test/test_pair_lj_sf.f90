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

integer(ib) :: Nsteps, Nprop, i
real(rb)    :: Rc, Rs, Rc2, Temp, Dt, Dt_2

type(tEmDee), target :: md
type(c_ptr), allocatable :: pair(:)

integer :: threads
character(256) :: filename, configFile

call command_line_arguments( filename, threads )
call read_data( filename )
call read_configuration( configFile )

allocate( pair(ntypes) )
do i = 1, ntypes
  pair(i) = EmDee_pair_lj_sf( epsilon(i), sigma(i) )
end do

call initialize_system( 1, pair )
call run( 0, Nprop )

contains
  include "common/contained.f90"
end program test

