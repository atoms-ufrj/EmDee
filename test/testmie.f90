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

program testfortran

use EmDee
use mConfig

implicit none

integer(ib) :: Nsteps, Nprop
real(rb)    :: Rc, Rs, Rc2, Temp, Dt, Dt_2
real(rb), pointer :: V(:,:)

integer(ib) :: step
type(tEmDee), target :: md
type(c_ptr), target :: pair, bond

integer :: i, j, argcount, threads
character(256) :: line, configFile

argcount = command_argument_count()
if (argcount == 1) then
  threads = 1
  call get_command_argument( 1, line )
else if (argcount == 2) then
  call get_command_argument( 1, line )
  read(line,*) threads
  call get_command_argument( 2, line )
else
  write(0,'("Usage: testfortran [number-of-threads] input-file")')
  stop
end if

call read_data( file = line )
call read_configuration( configFile )
call set_velocities

if (mod(N,2) /= 0) stop "PLEASE ENTER AN EVEN NUMBER OF ATOMS"
atomType(1:N/2) = 1
atomType(N/2+1:N) = 2
where (atomType == 1)
  Q = 1.0_rb
elsewhere
  Q = -1.0_rb
end where
Q(1:N/5) = 0.0_rb
Q(N/2+1:N/2+N/5) = 0.0_rb

!md = EmDee_system( threads, 1, Rc, Rs, N, c_loc(atomType), c_null_ptr )
md = EmDee_system( threads, 1, Rc, Rs, N, c_loc(atomType), c_null_ptr, c_null_ptr )

pair = EmDee_pair_mie_cut(1.0_rb, 1.0_rb,12.0_rb,6.0_rb)
!pair = EmDee_pair_lj_sf( 1.0_rb, 1.0_rb )
!pair = EmDee_pair_softcore_cut( 1.0_rb, 1.0_rb, 1.0_rb )

call EmDee_set_pair_model( md, 1, 1, pair, 1.0_rb )
call EmDee_set_pair_model( md, 2, 2, pair, 1.0_rb )
call EmDee_set_pair_model( md, 1, 2, pair, 1.0_rb )

call EmDee_set_coul_model( md, EmDee_coul_sf() )

call EmDee_upload( md, "charges"//c_null_char, c_loc(Q) )

do i = 1, N-1
  do j = i+1, N
    if (abs(i-j) < 10) call EmDee_ignore_pair( md, i, j )
  end do
end do

bond = EmDee_bond_harmonic( 1.0_rb, 1.0_rb )
call EmDee_add_bond( md, 1, 2, bond )
call EmDee_add_bond( md, 2, 3, bond )
call EmDee_add_bond( md, 4, 5, bond )

call EmDee_upload( md, "box"//c_null_char, c_loc(L) )
call EmDee_upload( md, "coordinates"//c_null_char, c_loc(R(1,1)) )
call EmDee_upload( md, "momenta"//c_null_char, c_loc(V(1,1)) )

print*, 0, md%Potential, md%Virial, md%Potential + md%Kinetic
do step = 1, Nsteps
  md%options%computeProps = mod(step,Nprop) == 0
  call EmDee_boost( md, 1.0_rb, 0.0_rb, Dt_2 )
  call EmDee_move( md, 1.0_rb, 0.0_rb, Dt )
  call EmDee_boost( md, 1.0_rb, 0.0_rb, Dt_2 )
  if (mod(step,Nprop) == 0) print*, step, md%Potential, md%Virial, md%Potential + md%Kinetic
end do
print*, "neighbor list builds = ", md%builds
print*, "pair time = ", md%pairTime, " s."
print*, "execution time = ", md%totalTime, " s."

contains
!---------------------------------------------------------------------------------------------------
  subroutine read_data( file )
    character(*), intent(in) :: file
    integer  :: inp, i, nseeds, seed
    open( newunit = inp, file = file, status = "old" )
    read(inp,*); read(inp,*) configFile
    read(inp,*); read(inp,*) Rc
    read(inp,*); read(inp,*) Rs
    read(inp,*); read(inp,*) seed
    read(inp,*); read(inp,*) Dt
    read(inp,*); read(inp,*) Nsteps
    read(inp,*); read(inp,*) Nprop
    read(inp,*); read(inp,*) Temp
    close(inp)
    call random_seed( size = nseeds )
    call random_seed( put = seed + 37*[(i-1,i=1,nseeds)] )
    Rc2 = Rc**2
    Dt_2 = 0.5_8*Dt
  end subroutine read_data
!---------------------------------------------------------------------------------------------------
  real(rb) function random_normal()
    real(rb) :: uni(2)
    call random_number( uni )
    random_normal = sqrt(-2.0_rb*log(uni(1))) * cos(6.283185307179586_rb*uni(2))
  end function random_normal
!---------------------------------------------------------------------------------------------------
subroutine create_configuration
  integer :: Nd, m, ind, i, j, k
  real(8) :: Vcm(3)
  Nd = ceiling(real(N,8)**(1.0_8/3.0_8))
  do ind = 1, N
    m = ind - 1
    k = m/(Nd*Nd)
    j = (m - k*Nd*Nd)/Nd
    i = m - j*Nd - k*Nd*Nd
    R(:,ind) = (L/real(Nd,8))*(real([i,j,k],8) + 0.5_8)
    V(:,ind) = [random_normal(), random_normal(), random_normal()]
  end do
  Vcm = sum(V,2)/N
  forall (i=1:N) V(:,i) = V(:,i) - Vcm
  V = sqrt(Temp*(3*N-3)/sum(V*V))*V
  Vcm = sum(V,2)/N
end subroutine create_configuration
!---------------------------------------------------------------------------------------------------
subroutine set_velocities
  integer :: i
  real(rb) :: Vcm(3)
  allocate( V(3,N) )
  do i = 1, N
    V(:,i) = [random_normal(), random_normal(), random_normal()]
  end do
  Vcm = sum(V,2)/N
  forall (i=1:N) V(:,i) = V(:,i) - Vcm
  V = sqrt(Temp*(3*N-3)/sum(V*V))*V
  Vcm = sum(V,2)/N
end subroutine set_velocities
!---------------------------------------------------------------------------------------------------
end program testfortran

