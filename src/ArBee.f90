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

module ArBee

use global
use math

implicit none

integer, parameter, private :: extra = 100

type tBody
  integer  :: NP        ! Number of particles
  integer  :: dof = 6   ! Number of degrees of freedom
  real(rb) :: mass      ! Total body mass
  real(rb) :: MoI(3)    ! Principal moments of inertia
  real(rb) :: rcm(3)    ! Center-of-mass position
  real(rb) :: pcm(3)    ! Center-of-mass momentum vector
  real(rb) :: q(0:3)    ! Unit quaternion of orientation
  real(rb) :: pi(0:3)   ! Quaternion momentum

  real(rb) :: omega(3)  ! angular velocities
  real(rb) :: F(3)      ! Resultant force
  real(rb) :: tau(3)    ! Resultant torque

  integer,  allocatable :: index(:)
  real(rb), allocatable :: M(:)
  real(rb), allocatable :: d(:,:)
  real(rb), allocatable :: delta(:,:)

  real(rb), allocatable :: E(:,:)
  real(rb), allocatable :: W(:,:)
  real(rb), allocatable :: F_R(:,:)

  real(rb) :: invMass    ! Inverse of body mass
  real(rb) :: invMoI(3)  ! Inverses of principal moments of inertia

  real(rb) :: I113       ! 1/(I1*(I1 - I3))
  real(rb) :: I313       ! 1/(I3*(I1 - I3))
  real(rb) :: I223       ! 1/(I2*(I2 - I3))
  real(rb) :: I212       ! 1/(I2*(I1 - I2))
  real(rb) :: m312       ! (I3 - I1)/I2

  contains

    procedure :: setup => tBody_setup
    procedure :: update => tBody_update
    procedure :: particle_momenta => tBody_particle_momenta
    procedure :: rotate_no_squish => tBody_rotate_no_squish
    procedure :: rotate_exact => tBody_rotate_exact
    procedure :: rotate_uniaxial => tBody_rotate_uniaxial
    procedure :: force_and_torque => tBody_force_and_torque
    procedure :: assign_momenta => tBody_assign_momenta

end type tBody

contains

!---------------------------------------------------------------------------------------------------

  pure subroutine tBody_setup( b, indexes, masses )
    class(tBody), intent(inout) :: b
    integer(ib),  intent(in)    :: indexes(:)
    real(rb),     intent(in)    :: masses(size(indexes))

    integer :: n

    n = size(indexes)
    b%NP = n
    allocate( b%index(n), b%M(n), b%d(3,n), b%delta(3,n), b%E(n,n), b%W(n,n), b%F_R(n,n) )
    b%index = indexes
    b%M = masses
    b%mass = sum(b%M)
    b%invMass = one/b%mass

  end subroutine tBody_setup

!---------------------------------------------------------------------------------------------------

  pure subroutine tBody_update( b, coords )
    class(tBody), intent(inout) :: b
    real(rb),     intent(in)    :: coords(3,b%NP)

    integer  :: x
    real(rb) :: inertia(3,3), A(3,3)

    ! Compute center-of-mass position and space-frame particle coordinates:
    forall (x=1:3) b%rcm(x) = sum(b%M*coords(x,:))*b%invMass
    forall (x=1:3) b%delta(x,:) = coords(x,:) - b%rcm(x)

    ! Compute upper-triangular inertia tensor:
    inertia(1,1) = sum(b%M*(b%delta(2,:)**2 + b%delta(3,:)**2))
    inertia(2,2) = sum(b%M*(b%delta(1,:)**2 + b%delta(3,:)**2))
    inertia(3,3) = sum(b%M*(b%delta(1,:)**2 + b%delta(2,:)**2))
    inertia(1,2) = -sum(b%M*b%delta(1,:)*b%delta(2,:))
    inertia(1,3) = -sum(b%M*b%delta(1,:)*b%delta(3,:))
    inertia(2,3) = -sum(b%M*b%delta(2,:)*b%delta(3,:))

    ! Diagonalize the inertia tensor:
    call diagonalization( inertia, A, b%MoI )
    A = transpose(A)
    b%invMoI = one/b%MoI
    associate(I1 => b%MoI(1), I2 => b%MoI(2), I3 => b%MoI(3))
      b%I113 = one/(I1*(I1 - I3))
      b%I313 = one/(I3*(I1 - I3))
      b%I223 = one/(I2*(I2 - I3))
      b%I212 = one/(I2*(I1 - I2))
      b%m312 = (I3 - I1)/I2
    end associate

    ! Compute quaternion:
    b%q = quaternion( A )

    ! Calculate position in the body-fixed frame:
    b%d = matmul( A, b%delta )

  end subroutine tBody_update

!---------------------------------------------------------------------------------------------------

  pure function matrix_B( q ) result( B )
    real(rb), intent(in) :: q(0:3)
    real(rb)             :: B(4,3)
    B = reshape( [-q(1),  q(0),  q(3), -q(2),  &
                  -q(2), -q(3),  q(0),  q(1),  &
                  -q(3),  q(2), -q(1),  q(0)], [4,3] )
  end function matrix_B

!---------------------------------------------------------------------------------------------------

  pure function matrix_C( q ) result( C )
    real(rb), intent(in) :: q(0:3)
    real(rb)             :: C(4,3)
    C = reshape( [-q(1),  q(0), -q(3),  q(2),  &
                  -q(2),  q(3),  q(0), -q(1),  &
                  -q(3), -q(2),  q(1),  q(0)], [4,3] )
  end function matrix_C

!---------------------------------------------------------------------------------------------------

  pure function matrix_Bt( q ) result( Bt )
    real(rb), intent(in) :: q(0:3)
    real(rb)             :: Bt(3,4)
    Bt = reshape( [-q(1), -q(2), -q(3), &
                    q(0), -q(3),  q(2), &
                    q(3),  q(0), -q(1), &
                   -q(2),  q(1),  q(0)  ], [3,4] )
  end function matrix_Bt

!---------------------------------------------------------------------------------------------------

  pure function matrix_Ct( q ) result( Ct )
    real(rb), intent(in) :: q(0:3)
    real(rb)             :: Ct(3,4)
    Ct = reshape( [-q(1), -q(2), -q(3), &
                    q(0),  q(3), -q(2), &
                   -q(3),  q(0),  q(1), &
                    q(2), -q(1),  q(0)  ], [3,4] )
  end function matrix_Ct

!---------------------------------------------------------------------------------------------------

  pure subroutine tBody_rotate_no_squish( b, delta_t, n )
    class(tBody), intent(inout) :: b
    real(rb),     intent(in)    :: delta_t
    integer(ib),  intent(in)    :: n
    integer :: i
    real(rb) :: dt, half_dt
    dt = delta_t/n
    half_dt = half*dt
    do i = 1, n
      call b % rotate_uniaxial( 3, half_dt )
      call b % rotate_uniaxial( 2, half_dt )
      call b % rotate_uniaxial( 1, dt )
      call b % rotate_uniaxial( 2, half_dt )
      call b % rotate_uniaxial( 3, half_dt )
    end do
    b%delta = matmul( matrix_Ct(b%q), matmul( matrix_B(b%q), b%d ) )
  end subroutine tBody_rotate_no_squish

!---------------------------------------------------------------------------------------------------

  pure subroutine tBody_rotate_uniaxial( b, k, dt )
    class(tBody), intent(inout) :: b
    integer,      intent(in)    :: k
    real(rb),     intent(in)    :: dt
    real(rb) :: BkQ(4), BkPi(4), omega_dt_by_2, vsin, vcos
    select case (k)
      case (1)
        BkQ  = [ -b%q(1),   b%q(0),   b%q(3),  -b%q(2)]
        BkPi = [-b%pi(1),  b%pi(0),  b%pi(3), -b%pi(2)]
      case (2)
        BkQ  = [ -b%q(2),  -b%q(3),   b%q(0),   b%q(1)]
        BkPi = [-b%pi(2), -b%pi(3),  b%pi(0),  b%pi(1)]
      case (3)
        BkQ  = [ -b%q(3),   b%q(2),  -b%q(1),   b%q(0)]
        BkPi = [-b%pi(3),  b%pi(2), -b%pi(1),  b%pi(0)]
    end select
    omega_dt_by_2 = dt*sum(b%pi*BkQ)/(4.0_rb*b%MoI(k))
    vsin = sin(omega_dt_by_2)
    vcos = cos(omega_dt_by_2)
    b%q  = vcos*b%q  + vsin*BkQ
    b%pi = vcos*b%pi + vsin*BkPi
  end subroutine tBody_rotate_uniaxial

!---------------------------------------------------------------------------------------------------

  pure subroutine tBody_rotate_exact( b, dt )

    class(tBody), intent(inout) :: b
    real(rb),     intent(in)    :: dt

    integer  :: i0, jump
    real(rb) :: w0(3), Iw(3), Lsq, TwoKr, r1, r3, l1, l3, lmin, a(3), m, K, wp, s0, c0, u0, u
    real(rb) :: L, deltaF, phi, z0(4), z(4), jac(3), inv2K

    w0 = b%omega
    Iw = b%MoI*w0
    Lsq = sum(Iw(2:3)**2)
    if (Lsq < epsilon(one)) then
      call b % rotate_uniaxial( 1, dt )
      return
    end if
    Lsq = Iw(1)**2 + Lsq
    L = sqrt(Lsq)
    TwoKr = sum(Iw*w0)
    r1 = Lsq - TwoKr*b%MoI(3)
    r3 = TwoKr*b%MoI(1) - Lsq
    l1 = b%I223*r1
    l3 = b%I212*r3
    lmin = min(l1,l3)
    a = [ sign(one,w0(1))*sqrt(b%I113*r1), sqrt(lmin), sign(one,w0(3))*sqrt(b%I313*r3) ]
    m = lmin/max(l1,l3)
    K = Carlson_RF( zero, one - m, one )
    inv2K = half/K
    s0 = w0(2)/a(2)
    if (abs(s0) < one) then
      c0 = merge( w0(1)/a(1), w0(3)/a(3), l1 < l3 )
      u0 = s0*Carlson_RF( one - s0*s0, one - m*s0*s0, one )
      i0 = staircase(u0*inv2K)
    else
      a(2) = abs(w0(2))
      s0 = sign(one,s0)
      c0 = zero
      u0 = sign(K,s0)
      i0 = 0
    end if
    wp = b%m312*a(1)*a(3)/a(2)
    u = wp*dt + u0
    jump = staircase(u*inv2K) - i0
    jac = jacobi( u, m )
    associate(sn => jac(1), cn => jac(2), dn => jac(3))
      if (l1 < l3) then
        b%omega = a*[cn,sn,dn]
        deltaF = deltaFcn( u0, c0, s0, w0(3)/a(3), u, cn, sn, dn, m, jump, b%MoI(1)*a(1)/L )
      else
        b%omega = a*[dn,sn,cn]
        deltaF = deltaFdn( u0, c0, s0, u, cn, sn, m, jump, b%MoI(1)*a(1)/L )
      end if
    end associate
    phi = (Lsq*(u - u0) + r3*deltaF)/(two*L*b%MoI(1)*wp)
    z0 = [Iw(3), Iw(2), L - Iw(1), zero]
    Iw = b%MoI*b%omega
    z = [ Iw(3), Iw(2), L - Iw(1), zero]*cos(phi) + &
        [-Iw(2), Iw(3), zero, L - Iw(1)]*sin(phi)
    b%q = normalize( z*sum(z0*b%q) + matmul( matrix_C(z), matmul( matrix_Ct(z0), b%q ) ) )
    b%pi = matmul( matrix_B(b%q), two*Iw )
    b%delta = matmul( matrix_Ct(b%q), matmul( matrix_B(b%q), b%d ) )
    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      pure function Theta( x, n, m )
        real(rb), intent(in) :: x, n, m
        real(rb)             :: Theta
        real(rb) :: x2
        x2 = x*x
        Theta = -third*n*x*x2*Carlson_RJ( one - x2, one - m*x2, one, one + n*x2 )
      end function Theta
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      pure function deltaFcn( u0, c0, s0, d0, u, cn, sn, dn, m, jump, alpha )
        integer,  intent(in) :: jump
        real(rb), intent(in) :: u0, c0, s0, d0, u, cn, sn, dn, m, alpha
        real(rb)             :: deltaFcn
        real(rb) :: eta, C
        eta = alpha**2/(one - alpha**2)
        C = sqrt(m + eta)
        deltaFcn = u - u0 + sign(one,cn)*Theta(sn,eta,m) - sign(one,c0)*Theta(s0,eta,m) &
                          + (alpha/C)*(atan(C*sn/dn) - atan(C*s0/d0))
        if (jump /= 0) deltaFcn = deltaFcn + jump*two*Theta(one,eta,m)
        deltaFcn = (eta + one)*deltaFcn
      end function deltaFcn
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      pure function deltaFdn( u0, c0, s0, u, cn, sn, m, jump, alpha )
        integer,  intent(in) :: jump
        real(rb), intent(in) :: u0, c0, s0, u, cn, sn, m, alpha
        real(rb)             :: deltaFdn
        real(rb) :: eta, k2eta, C
        eta = alpha**2
        eta = eta/(one - eta)
        k2eta = m*eta
        C = sqrt(one + k2eta)
        deltaFdn = u - u0 + sign(one,cn)*Theta(sn,k2eta,m) - sign(one,c0)*Theta(s0,k2eta,m) &
                          + (alpha/C)*(atan(C*sn/cn) - atan(C*s0/c0))
        if (jump /= 0) deltaFdn = deltaFdn + jump*(two*Theta(one,k2eta,m) + (alpha/C)*Pi)
        deltaFdn = (eta + one)*deltaFdn
      end function deltaFdn
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine tBody_rotate_exact

!---------------------------------------------------------------------------------------------------

  pure function tBody_particle_momenta( b ) result( P )
    class(tBody), intent(in) :: b
    real(rb)                 :: P(3,b%NP)
    integer  :: k
    real(rb) :: omega(3), At(3,3)
    At = matmul( matrix_Ct(b%q), matrix_B(b%q) )
    omega = matmul(At, b%omega)
    forall(k=1:b%NP) P(:,k) = b%M(k)*(b%invMass*b%pcm + cross_product(omega, b%delta(:,k)))
  end function tBody_particle_momenta

!---------------------------------------------------------------------------------------------------

  pure subroutine tBody_force_and_torque( b, F )
    class(tBody), intent(inout) :: b
    real(rb),     intent(in)    :: F(:,:)
    integer :: j
    real(rb) :: Fj(3)
    b%F = zero
    b%tau = zero
    do j = 1, b%NP
      Fj = F(:,b%index(j))
      b%F = b%F + Fj
      b%tau = b%tau + cross_product( b%delta(:,j), Fj )
    end do
  end subroutine tBody_force_and_torque

!---------------------------------------------------------------------------------------------------

  pure subroutine tBody_assign_momenta( b, input )
    class(tBody), intent(inout) :: b
    real(rb),     intent(in)    :: input(:)
    if (size(input) == 3) then
      b%omega = input
      b%pi = matmul( matrix_B(b%q), two*b%MoI*b%omega )
    else
      b%pi = input
      b%omega = half*b%invMoI*matmul( matrix_Bt(b%q), b%pi )
    end if
  end subroutine tBody_assign_momenta

!---------------------------------------------------------------------------------------------------

end module ArBee
