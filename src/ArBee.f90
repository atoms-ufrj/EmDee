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

  real(rb) :: invMass    ! Inverse of body mass
  real(rb) :: invMoI(3)  ! Inverses of principal moments of inertia

  contains

    procedure :: setup => tBody_setup
    procedure :: update => tBody_update
    procedure :: particle_momenta => tBody_particle_momenta
!    procedure :: rotate => tBody_rotate_no_squish
    procedure :: rotate => tBody_rotate_analytical
    procedure :: force_torque_virial => tBody_force_torque_virial
    procedure :: assign_momenta => tBody_assign_momenta

end type tBody

contains

!---------------------------------------------------------------------------------------------------

  subroutine tBody_setup( b, indexes, masses )
    class(tBody), intent(inout) :: b
    integer(ib),  intent(in)    :: indexes(:)
    real(rb),     intent(in)    :: masses(size(indexes))

    b%NP = size(indexes)
    allocate( b%index(b%NP), b%M(b%NP), b%d(3,b%NP), b%delta(3,b%NP) )
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

  elemental subroutine tBody_rotate_no_squish( b, dt )
    class(tBody), intent(inout) :: b
    real(rb),     intent(in)    :: dt
    real(rb) :: half_dt
    half_dt = half*dt
    call uniaxial_rotation( b, 3, half_dt )
    call uniaxial_rotation( b, 2, half_dt )
    call uniaxial_rotation( b, 1, dt )
    call uniaxial_rotation( b, 2, half_dt )
    call uniaxial_rotation( b, 3, half_dt )
    b%delta = matmul( matmul( matrix_Ct(b%q), matrix_B(b%q) ), b%d )
    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      elemental subroutine uniaxial_rotation( b, k, dt )
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
      end subroutine uniaxial_rotation
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine tBody_rotate_no_squish

!---------------------------------------------------------------------------------------------------

  subroutine tBody_rotate_analytical( b, dt )
    class(tBody), intent(inout) :: b
    real(rb),     intent(in)    :: dt

    integer  :: status
    real(rb) :: w0(3), Iw(3), Lsq, TwoKr, l1, l3, a(3), ksq, wp, s0, u0, u
    real(rb) :: L, deltaF, phi, z0(4), z(4)
    real(rb), target :: cn, sn, dn

    w0 = b%omega
    Iw = b%MoI*w0
    TwoKr = sum(Iw*w0)
    Lsq = sum(Iw*Iw)
    L = sqrt(Lsq)
    z0 = [Iw(3), Iw(2), L - Iw(1), zero]/sqrt(two*L*(L - Iw(1)))
    associate (I1 => b%MoI(1), I2 => b%MoI(2), I3 => b%MoI(3))
      l1 = sqrt((Lsq - TwoKr*I3)/(I2*(I2 - I3)))
      l3 = sqrt((TwoKr*I1 - Lsq)/(I2*(I1 - I2)))
      a = [ w0(1)/sqrt(one - (w0(2)/l1)**2), min(l1,l3), w0(3)/sqrt(one - (w0(2)/l3)**2) ]
      ksq = (a(2)/max(l1,l3))**2
      wp = (I3 - I1)*a(1)*a(3)/(I2*a(2))
      s0 = w0(2)/a(2)
      u0 = s0*gsl_sf_ellint_RF( one - s0*s0, one - ksq*s0*s0, one, GSL_PREC_DOUBLE )
      u = wp*dt + u0
      status = gsl_sf_elljac_e( u, ksq, c_loc(sn), c_loc(cn), c_loc(dn) )
      if (l1 < l3) then
        b%omega = a*[cn,sn,dn]
        deltaF = deltaFcn( u0, w0(1)/a(1), s0, w0(3)/a(3), u, cn, sn, dn, ksq, I1*a(1)/L )
      else
        b%omega = a*[dn,sn,cn]
        deltaF = deltaFdn( u0, w0(3)/a(3), s0, u, cn, sn, ksq, I1*a(1)/L )
      end if
      phi = (Lsq*(u - u0) + (TwoKr*I1 - Lsq)*deltaF)/(two*L*I1*wp)
    end associate
    Iw = b%MoI*b%omega
    z = ([ Iw(3), Iw(2), L - Iw(1), zero]*cos(phi) + &
         [-Iw(2), Iw(3), zero, L - Iw(1)]*sin(phi) )/sqrt(two*L*(L - Iw(1)))
    b%q = z*sum(z0*b%q) + matmul(matrix_C(z),matmul(matrix_Ct(z0),b%q))
    b%pi = matmul( matrix_B(b%q), two*b%MoI*b%omega )
    b%delta = matmul( matmul( matrix_Ct(b%q), matrix_B(b%q) ), b%d )

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      function Theta( x, n, m )
        real(rb), intent(in) :: x, n, m
        real(rb)             :: Theta
        real(rb) :: xsq
        xsq = x*x
        Theta = -third*n*x*xsq* &
                 gsl_sf_ellint_RJ( one - xsq, one - m*xsq, one, one + n*xsq, GSL_PREC_DOUBLE )
      end function Theta
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      function deltaFcn( u0, c0, s0, d0, u, cn, sn, dn, ksq, alpha )
        real(rb), intent(in) :: u0, c0, s0, d0, u, cn, sn, dn, ksq, alpha
        real(rb)             :: deltaFcn
        integer  :: jump
        real(rb) :: eta, C, inv2K
        eta = alpha**2/(one - alpha**2)
        C = sqrt(ksq + eta)
        deltaFcn = u - u0 + sign(one,cn)*Theta(sn,eta,ksq) - sign(one,c0)*Theta(s0,eta,ksq) &
                          + (alpha/C)*(atan(C*sn/dn) - atan(C*s0/d0))
        inv2K = half/gsl_sf_ellint_RF( zero, one - ksq, one, GSL_PREC_DOUBLE )
        jump = nint(u*inv2K) - nint(u0*inv2K)
        if (jump /= 0) deltaFcn = deltaFcn + jump*two*Theta(one,eta,ksq)
        deltaFcn = (eta + one)*deltaFcn
      end function deltaFcn
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      function deltaFdn( u0, c0, s0, u, cn, sn, ksq, alpha )
        real(rb), intent(in) :: u0, c0, s0, u, cn, sn, ksq, alpha
        real(rb)             :: deltaFdn
        integer  :: jump
        real(rb) :: eta, k2eta, C, inv2K
        eta = alpha**2/(one - alpha**2)
        k2eta = ksq*eta
        C = sqrt(one + k2eta)
        deltaFdn = u - u0 + sign(one,cn)*Theta(sn,k2eta,ksq) - sign(one,c0)*Theta(s0,k2eta,ksq) &
                          + (alpha/C)*(atan(C*sn/cn) - atan(C*s0/c0))
        inv2K = half/gsl_sf_ellint_RF( zero, one - ksq, one, GSL_PREC_DOUBLE )
        jump = nint(u*inv2K) - nint(u0*inv2K)
        if (jump /= 0) deltaFdn = deltaFdn + jump*(two*Theta(one,k2eta,ksq) + (alpha/C)*Pi)
        deltaFdn = (eta + one)*deltaFdn
      end function deltaFdn
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine tBody_rotate_analytical

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

  function tBody_force_torque_virial( b, F ) result( virial )
    class(tBody), intent(inout) :: b
    real(rb),     intent(in)    :: F(:,:)
    real(rb)                    :: virial
    integer :: j
    real(rb) :: Fj(3)
    b%F = zero
    b%tau = zero
    virial = zero
    do j = 1, b%NP
      Fj = F(:,b%index(j))
      b%F = b%F + Fj
      b%tau = b%tau + cross_product( b%delta(:,j), Fj )
      virial = virial + sum(b%delta(:,j)*Fj)
    end do
  end function tBody_force_torque_virial

!---------------------------------------------------------------------------------------------------

  subroutine tBody_assign_momenta( b, input )
    class(tBody), intent(inout) :: b
    real(rb),     intent(in)    :: input(:)
    if (size(input) == 3) then
      b%omega = input
      b%pi = matmul( matrix_B(b%q), two*b%MoI*b%omega )
    else if (size(input) == 4) then
      b%pi = input
      b%omega = half*b%invMoI*matmul( matrix_Bt(b%q), b%pi )
    else
      stop "ERROR: invalid rigid body momentum assignment"
    end if
  end subroutine tBody_assign_momenta

!---------------------------------------------------------------------------------------------------

end module ArBee
