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

program testfortran

use EmDee

implicit none

integer(ib) :: N, Nsteps, Nprop
real(rb)    :: rho, L, Rc, Rs, Rc2, Temp, Dt, Dt_2
real(rb), pointer :: R(:,:), V(:,:), F(:,:)

integer(ib) :: step
real(rb)    :: ti, tf
type(c_ptr) :: mdp
type(tEmDee), target :: md
type(tModel), target :: lj, bond

integer :: i, j
!integer, pointer :: first(:), last(:), item(:)

call read_data
call create_configuration
mdp = c_loc(md)
call md_initialize( mdp, Rc, Rs, N, 1, c_null_ptr, c_null_ptr, c_loc(R(1,1)), c_loc(F(1,1)) )

lj = pair_lj( 1.0_rb, 1.0_rb )
call md_set_pair( c_loc(md), 1, 1, c_loc(lj) )

!print*, match( [1,7,5,3,2,4,8,6], [2,3,4,8])
!stop

do i = 1, N-1
  do j = i+1, N
    if (abs(i-j) < 10) call md_exclude_pair( mdp, i, j )
  end do
end do

bond = bond_harmonic( 1.0_rb, 1.0_rb )
call md_add_bond( mdp, 1, 2, c_loc(bond) )
call md_add_bond( mdp, 2, 3, c_loc(bond) )
call md_add_bond( mdp, 4, 5, c_loc(bond) )

!call c_f_pointer( md%excluded%first, first, [md%natoms])
!call c_f_pointer( md%excluded%last, last, [md%natoms])
!call c_f_pointer( md%excluded%item, item, [md%excluded%nitems])

!do i = 1, N
!  print*, i, " --- ", item(first(i):last(i))
!end do
!stop


call md_compute_forces( mdp, L )
print*, 0, md%Energy, md%Virial
call cpu_time( ti )
do step = 1, Nsteps
  if (mod(step,Nprop) == 0) print*, step, md%Energy, md%Virial
  V = V + Dt_2*F
  R = R + Dt*V
  call md_compute_forces( mdp, L )
  V = V + Dt_2*F
end do
call cpu_time( tf )
print*, "neighbor list builds = ", md%builds
print*, "pair time = ", md%time, " s."
print*, "execution time = ", tf - ti, " s."

contains
!---------------------------------------------------------------------------------------------------
  subroutine read_data
    integer  :: i, nseeds, seed
    read(*,*); read(*,*) N
    read(*,*); read(*,*) Rc
    read(*,*); read(*,*) Rs
    read(*,*); read(*,*) seed
    read(*,*); read(*,*) Dt
    read(*,*); read(*,*) Nsteps
    read(*,*); read(*,*) Nprop
    read(*,*); read(*,*) rho
    read(*,*); read(*,*) Temp
    call random_seed( size = nseeds )
    call random_seed( put = seed + 37*[(i-1,i=1,nseeds)] )
    Rc2 = Rc**2
    L = (N/rho)**(1.0_8/3.0_8)
    Dt_2 = 0.5_8*Dt
    allocate( R(3,N), V(3,N), F(3,N) )
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
end program testfortran

