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

md = EmDee_system( threads, 1, Rc, Rs, N, c_loc(atomType), c_loc(mass), c_loc(molecule) )
do i = 1, ntypes
  if (epsilon(i) == 0.0_rb) then
    pair = EmDee_pair_none( )
  else
    pair = EmDee_pair_lj_cut( epsilon(i), sigma(i) )
  end if
  call EmDee_set_pair_model( md, i, i, pair, kCoul )
end do
call EmDee_set_coul_model( md, EmDee_coul_long() )
!call EmDee_set_coul_model( md, EmDee_coul_none() )
call EmDee_set_kspace_model( md, EmDee_kspace_ewald( 0.0001_rb ) )
call EmDee_upload( md, "charges"//c_null_char, c_loc(Q) )
call EmDee_upload( md, "coordinates"//c_null_char, c_loc(R(1,1)) )
call EmDee_upload( md, "box"//c_null_char, c_loc(L) )
call EmDee_random_momenta( md, kB*Temp, .true._1, seed )

md%Options%RotationMode = 1
call run( Nsteps, Nprop )
call runxyz( 5000, 10 )

contains
  include "common/contained.f90"
!---------------------------------------------------------------------------------------------------
  subroutine runxyz( Nsteps, Nprop )
    integer, intent(in) :: Nsteps, Nprop

    integer :: step

    open(unit=10,file="rods.xyz",status="replace")
    do step = 1, Nsteps
      md%Energy%Compute = mod(step,Nprop) == 0
      call EmDee_advance( md, 0.0_rb, 0.0_rb, dt )
      if (mod(step,Nprop) == 0) then
        call EmDee_download( md, "coordinates"//char(0), c_loc(R(1,1)) )
        write(10,*) N
        write(10,*)
        do i = 1, N
          write(10,*) atomType(i), R(:,i)
        end do
      end if
    end do
    close(10)

  end subroutine runxyz
!---------------------------------------------------------------------------------------------------
end program test
