program testfortran

use mEmDee
use mRandom
use, intrinsic :: iso_c_binding

implicit none

integer, parameter :: ib = c_int, &
                      rb = c_double

integer(ib) :: N, seed, passo, Npassos, Nprop
real(rb)    :: rho, L, Rc, Rs, Rc2, Temp, Press, Dt, InvL, Ws, SixWs, Ec, Dt_2
real(rb), pointer :: R(:,:), V(:,:), F(:,:)
real(rb), target :: mass
type(shr3) :: random
type(tEmDee), target :: md

integer :: i
real(8) :: ti, tf
type(c_ptr) :: mdp

call read_data
call create_configuration
mdp = c_loc(md)
mass = 1.0_rb
call md_initialize( mdp, Rc, Rs, N, 1, c_null_ptr, c_loc(mass) )
call md_set_lj( mdp, 1, 1, 1.0_rb, 1.0_rb )
call md_upload( mdp, c_loc(R), c_loc(V) )
call md_compute_forces( mdp, L )
print*, md%Energy, md%Virial
print*, md%npairs
print*, "L = ", L
call cpu_time( ti )
do passo = 1, Npassos
  if (mod(passo,Nprop) == 0) print*, passo, md%Energy
  call md_change_momenta( mdp, 1.0_rb, Dt_2 )
  call md_change_coordinates( mdp, 1.0_rb, Dt )
  call md_compute_forces( mdp, L )
  call md_change_momenta( mdp, 1.0_rb, Dt_2 )
end do
call cpu_time( tf )
print*, "npairs = ", md%npairs
print*, "execution time = ", tf - ti, " s."
contains
!---------------------------------------------------------------------------------------------------
subroutine read_data
  integer  :: nchain, nloop, nsy
  real(8) :: InvRc2, InvRc6, InvRc12, tdamp
  read(*,*); read(*,*) N
  read(*,*); read(*,*) rho
  read(*,*); read(*,*) Rc, Rs
  read(*,*); read(*,*) Temp
  read(*,*); read(*,*) Press
  read(*,*); read(*,*) seed
  read(*,*); read(*,*) Dt
  read(*,*); read(*,*) Npassos
  read(*,*); read(*,*) Nprop
  read(*,*); read(*,*) nchain, nloop, nsy
  read(*,*); read(*,*) tdamp
  call random % setup( seed )
  Rc2 = Rc**2
  L = (N/rho)**(1.0_8/3.0_8)
  InvL = 1.0_8/L
  InvRc2 = 1.0_8/Rc2
  InvRc6 = InvRc2*InvRc2*InvRc2
  InvRc12 = InvRc6*InvRc6
  Ws = 2.0_8*InvRc12 - InvRc6
  SixWs = 6.0_8*Ws
  Ec = InvRc12 - InvRc6
  Dt_2 = 0.5_8*Dt
  allocate( R(3,N), V(3,N), F(3,N) )
end subroutine read_data
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
    V(:,ind) = [random%normal(), random%normal(), random%normal()]
  end do
  Vcm = sum(V,2)/N
  forall (i=1:N) V(:,i) = V(:,i) - Vcm
  V = sqrt(Temp*(3*N-3)/sum(V*V))*V
  Vcm = sum(V,2)/N
end subroutine create_configuration
!---------------------------------------------------------------------------------------------------
end program testfortran

