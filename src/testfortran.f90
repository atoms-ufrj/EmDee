program testfortran

use mEmDee
use mRandom
use, intrinsic :: iso_c_binding

implicit none

integer, parameter :: ib = c_int, &
                      rb = c_double

integer(ib) :: N, seed, Npassos, Nprop
real(rb)    :: rho, L, Rc, Rs, Rc2, Temp, Press, Dt, InvL, Ws, SixWs, Ec, Dt_2
real(rb), pointer :: R(:,:), V(:,:), F(:,:)
type(shr3) :: random
type(tEmDee), target :: md

integer :: i
real(8) :: ti, tf
type(c_ptr) :: mdp

!integer :: icell, ix, iy, iz, m, mm, count
!character(3) :: adv

!n = 3
!m = 2*n + 1
!mm = m*m
!count = 0
!do icell = (m*mm-1)/2+1, m*mm-1
!  iz = icell/mm
!  iy = (icell - iz*mm)/m
!  ix = icell - iy*m - iz*mm
!  Dt = sqrt(real(max(0,abs(ix-n)-1)**2 + max(0,abs(iy-n)-1)**2 + max(0,abs(iz-n)-1)**2,rb))
!  if (Dt < real(n,rb)) then
!    count = count + 1
!    if (mod(count,5) == 0) then
!      adv = "yes"
!    else
!      adv = "no"
!    end if
!    write(*,'(" {",I2,",",I2,",",I2,"},")',advance=adv) ix-n, iy-n, iz-n
!  end if
!!  print*, icell, Dt
!end do



!do k = -n, 0
!  do j = -n, 0
!    do i = -n, 0
!      Dt = sqrt(real(max(0,abs(i)-1)**2 + max(0,abs(j)-1)**2 + max(0,abs(k)-1)**2,rb))
!      if (Dt < real(n,rb)) count = count + 1
!      print*, (i+n)+m*((j+n)+m*(k+n)), Dt, Dt < real(n,rb)
!    end do
!  end do
!end do
!print*
!print*, count
!stop

call read_data
call create_configuration
mdp = c_loc(md)
call md_initialize( mdp, Rc, Rs, N, c_null_ptr )
call md_set_lj( mdp, 1, 1, 1.0_rb, 1.0_rb )
call md_upload( mdp, c_loc(R), c_loc(V) )
call md_compute_forces( mdp, L )
print*, md%npairs
print*, md%Energy
call cpu_time( ti )
do i = 1, 200
  call md_compute_forces( mdp, L )
end do
!do i = 1, 5000
!  print*, i, md%Energy
!  call md_change_momenta( mdp, 1.0_rb, Dt_2 )
!  call md_change_coordinates( mdp, 1.0_rb, Dt )
!  call md_compute_forces( mdp, L )
!  call md_change_momenta( mdp, 1.0_rb, Dt_2 )
!end do
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

