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
  pair = EmDee_shifted_force( EmDee_pair_lj_cut( epsilon(i), sigma(i) ) )
  call EmDee_set_pair_model( md, i, i, pair, kCoul )
end do
call EmDee_upload( md, "charges"//c_null_char, c_loc(Q) )
call EmDee_upload( md, "box"//c_null_char, c_loc(L) )
call EmDee_upload( md, "coordinates"//c_null_char, c_loc(R(1,1)) )
call set_momenta
call EmDee_upload( md, "momenta"//c_null_char, c_loc(P(1,1)) )
call run( Nsteps, Nprop )

contains
!---------------------------------------------------------------------------------------------------
  real(rb) function random_normal()
    real(rb) :: uni(2)
    call random_number( uni )
    random_normal = sqrt(-2.0_rb*log(uni(1))) * cos(6.283185307179586_rb*uni(2))
  end function random_normal
!---------------------------------------------------------------------------------------------------
  subroutine set_momenta
    integer :: i
    real(rb) :: Vcm(3), V(3,N), M(N)
    do i = 1, N
      V(:,i) = [random_normal(), random_normal(), random_normal()]
    end do
    M = mass(atomType)
    forall (i=1:3) Vcm(i) = sum(M*V(i,:))/sum(M)
    forall (i=1:N) V(:,i) = V(:,i) - Vcm
    V = sqrt((3*N-3)*kB*Temp/sum(M*sum(V*V,1)))*V
    forall (i=1:N) P(:,i) = M(i)*V(:,i)
  end subroutine set_momenta
!---------------------------------------------------------------------------------------------------
  include "common/contained.f90"
end program test
