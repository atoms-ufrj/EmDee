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

use global
use kspace_ewald_module

implicit none

#include "emdee.f03"

real(rb), parameter :: mvv2e = 2390.057364_rb         ! Da*A²/fs² to kcal/mol
real(rb), parameter :: kB = 8.3144621455E-7_rb        ! Boltzmann constant in Da*A²/(fs²*K)
real(rb), parameter :: Pconv = 1.6388244954E+8_rb     ! Da/(A*fs²) to atm
real(rb), parameter :: kCoul = 0.13893545755135628_rb ! Coulomb constant in Da*A³/(fs*e)²

character(sl) :: configFile
integer(ib)   :: ntypes, NperMol, N, Nmol
integer(ib), pointer :: atomType(:)
real(rb),    pointer :: atomMass(:), Q(:)
real(rb),    allocatable :: mass(:), eps(:), sigma(:), charge(:)

integer(ib) :: i, Nsteps, Nprop
real(rb)    :: rho, Rc, Rs, Rc2, Temp, Dt, Dt_2, Elong
real(rb), target  :: L(3)
real(rb), pointer :: R(:,:)

type(tEmDee), target  :: md
type(c_ptr),  pointer :: pair(:)

type(kspace_ewald) :: kspace
!type(c_ptr) :: bond

integer :: argcount, threads
character(256) :: line

argcount = command_argument_count()
if (argcount == 1) then
  threads = 1
  call get_command_argument( 1, line )
else if (argcount == 2) then
  call get_command_argument( 1, line )
  read(line,*) threads
  call get_command_argument( 2, line )
else
  write(0,'("Usage: test_kspace_ewald [number-of-threads] input-file")')
  stop
end if

call read_data( file = line )
call read_configuration

md = EmDee_system( threads, 1, Rc, Rs, N, c_loc(atomType), c_loc(atomMass) )

allocate( pair(ntypes) )
do i = 1, ntypes
  if (eps(i) == zero) then
    pair(i) = EmDee_pair_none()
  else
    pair(i) = EmDee_pair_lj_cut( eps(i), sigma(i) )
  end if
  call EmDee_set_pair_model( md, i, i, pair(i) )
end do

!bond = EmDee_bond_none( )
do i = 1, Nmol
!  call EmDee_add_bond( md, (i-1)*NperMol+1, (i-1)*NperMol+2, bond )
!  call EmDee_add_bond( md, (i-1)*NperMol+2, (i-1)*NperMol+3, bond )
  call EmDee_ignore_pair( md, (i-1)*NperMol+1, (i-1)*NperMol+2 )
  call EmDee_ignore_pair( md, (i-1)*NperMol+2, (i-1)*NperMol+3 )
end do

call EmDee_upload( md, "charges"//c_null_char, c_loc(Q) )
call EmDee_upload( md, "box"//c_null_char, c_loc(L) )
call EmDee_upload( md, "coordinates"//c_null_char, c_loc(R(1,1)) )
print*, md%Potential/kB, md%Virial/kB

call kspace%setup( [5.6_rb], [5] )
call kspace%initialize( 1, Rc, Q )
call kspace%compute( 1, Q, R/L(1), L(1)**3, Elong )
print*, Elong/kB

contains
!---------------------------------------------------------------------------------------------------
  subroutine read_data( file )
    character(*), intent(in) :: file
    integer  :: inp, i, nseeds, seed
    open( newunit = inp, file = file, status = "old" )
    read(inp,*); read(inp,*) configFile
    read(inp,*); read(inp,*) ntypes
    allocate( mass(ntypes), eps(ntypes), sigma(ntypes), charge(ntypes) )
    read(inp,*)
    do i = 1, ntypes
      read(inp,*) mass(i), eps(i), sigma(i), charge(i)
    end do
    eps = kB*eps
    read(inp,*); read(inp,*) Rc
    read(inp,*); read(inp,*) Rs
    read(inp,*); read(inp,*) seed
    read(inp,*); read(inp,*) Dt
    read(inp,*); read(inp,*) Nsteps
    read(inp,*); read(inp,*) Nprop
    read(inp,*); read(inp,*) rho
    read(inp,*); read(inp,*) Temp
    close(inp)
    call random_seed( size = nseeds )
    call random_seed( put = seed + 37*[(i-1,i=1,nseeds)] )
    Rc2 = Rc**2
    Dt_2 = 0.5_8*Dt
    allocate( R(3,N) )
  end subroutine read_data
!---------------------------------------------------------------------------------------------------
  real(rb) function random_normal()
    real(rb) :: uni(2)
    call random_number( uni )
    random_normal = sqrt(-2.0_rb*log(uni(1)))*cos(6.283185307179586_rb*uni(2))
  end function random_normal
!---------------------------------------------------------------------------------------------------
subroutine read_configuration
  integer :: inp, i, id
  open( newunit = inp, file = configFile, status = "old" )
  read(inp,*); read(inp,*) L
  read(inp,*); read(inp,*) NperMol
  read(inp,*); read(inp,*) Nmol
  N = NperMol*Nmol
  allocate( R(3,N), atomType(N), atomMass(N), Q(N) )
  read(inp,*)
  do i = 1, N
    read(inp,*) id, R(:,id), atomType(id)
  end do
  close(inp)
  atomMass = mass(atomType)
  Q = charge(atomType)
end subroutine read_configuration
!---------------------------------------------------------------------------------------------------
end program test

