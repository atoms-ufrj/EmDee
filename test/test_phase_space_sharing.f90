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

type(tEmDee) :: other
integer :: step

call command_line_arguments( filename, threads )
call read_data( filename )
call read_configuration( configFile )
call unit_conversions

md = EmDee_system( threads, 1, Rc, Rs, N, c_loc(atomType), c_loc(mass), c_loc(molecule) )
other = EmDee_system( threads, 1, 0.5*Rc, Rs, N, c_loc(atomType), c_loc(mass), c_loc(molecule) )
do i = 1, ntypes
  if (epsilon(i) == 0.0_rb) then
    pair = EmDee_pair_none( )
  else
    pair = EmDee_shifted_force( EmDee_pair_lj_cut( epsilon(i), sigma(i) ) )
  end if
  call EmDee_set_pair_model( md, i, i, pair, kCoul )
  call EmDee_set_pair_model( other, i, i, pair, kCoul )
end do

call EmDee_set_coul_model( md, EmDee_coul_sf() )
call EmDee_upload( md, "charges"//c_null_char, c_loc(Q) )
call EmDee_upload( md, "coordinates"//c_null_char, c_loc(R(1,1)) )
call EmDee_upload( md, "box"//c_null_char, c_loc(L) )

call EmDee_set_coul_model( other, EmDee_coul_sf() )
call EmDee_upload( other, "charges"//c_null_char, c_loc(Q) )
call EmDee_upload( other, "coordinates"//c_null_char, c_loc(R(1,1)) )
call EmDee_upload( other, "box"//c_null_char, c_loc(L) )

call run( 100, Nprop )

call EmDee_share_phase_space( md, other )
other%Options%Compute = .true.
call EmDee_compute_forces( other )

print*
associate( system => other )
  print*, 0, mvv2e*system%Energy%Potential, mvv2e*system%Virial%Total, &
        mvv2e*(system%Energy%Potential + system%Kinetic%Total)
  do step = 1, 100
    system%Options%Compute = mod(step,Nprop) == 0
    call EmDee_boost( md, 1.0_rb, 0.0_rb, Dt_2 )
    call EmDee_displace( md, 1.0_rb, 0.0_rb, Dt )
    call EmDee_boost( md, 1.0_rb, 0.0_rb, Dt_2 )
    if (mod(step,Nprop) == 0) then
      print*, step, mvv2e*system%Energy%Potential, mvv2e*system%Virial%Total, &
              mvv2e*(system%Energy%Potential + system%Kinetic%Total)
    end if
  end do
  print*, "neighbor list builds = ", system%builds
  print*, "pair time      = ", system%Time%Pair, " s."
  print*, "motion time    = ", system%Time%Motion, " s."
  print*, "neighbor time  = ", system%Time%Neighbor, " s."
  print*, "execution time = ", system%Time%Total, " s."
end associate

contains
  include "common/contained.f90"
end program test
