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

type rigidBody
  integer  :: NP        ! Number of particles
  integer  :: dof = 6   ! Number of degrees of freedom
  real(rb) :: mass      ! Total body mass
  real(rb) :: MoI(3)    ! Principal moments of inertia
  real(rb) :: rcm(3)    ! Center-of-mass position
  real(rb) :: pcm(3)    ! Center-of-mass momentum vector
  real(rb) :: q(0:3)    ! Unit quaternion of orientation
  real(rb) :: pi(0:3)   ! Quaternion momentum

  real(rb) :: F(3)      ! Resultant force
  real(rb) :: tau(3)    ! Resultant torque

  integer,  allocatable :: index(:)
  real(rb), allocatable :: M(:)
  real(rb), allocatable :: d(:,:)
  real(rb), allocatable :: delta(:,:)

  real(rb) :: invMass    ! Inverse of body mass
  real(rb) :: invMoI(3)  ! Inverses of principal moments of inertia

  contains

    procedure :: setup => rigidBody_setup
    procedure :: update => rigidBody_update
    procedure :: particle_momenta => rigidBody_particle_momenta
    procedure :: rotate => rigidBody_rotate
    procedure :: force_torque_virial => rigidBody_force_torque_virial

end type rigidBody

contains

!---------------------------------------------------------------------------------------------------

  subroutine realloc_rigid_body_list( list, Nmax )
    type(c_ptr), intent(inout) :: list
    integer(ib), intent(inout) :: Nmax

    type(rigidBody), pointer :: old(:), new(:)

    allocate( new(Nmax+extra) )
    if (c_associated(list)) then
      call c_f_pointer( list, old, [Nmax] )
      new(1:Nmax) = old
      deallocate( old )
    end if
    Nmax = Nmax + extra
    list = c_loc(new)

  end subroutine realloc_rigid_body_list

!---------------------------------------------------------------------------------------------------

  subroutine rigidBody_setup( b, indexes, masses )
    class(rigidBody), intent(inout) :: b
    integer(ib),      intent(in)    :: indexes(:)
    real(rb),         intent(in)    :: masses(size(indexes))

    b%NP = size(indexes)
    allocate( b%index(b%NP), b%M(b%NP), b%d(3,b%NP), b%delta(3,b%NP) )
    b%index = indexes
    b%M = masses
    b%mass = sum(b%M)
    b%invMass = one/b%mass
    b%rcm = zero
    b%d = zero
    b%delta = zero
    b%q = zero
    b%pcm = zero
    b%pi  = zero

  end subroutine rigidBody_setup

!---------------------------------------------------------------------------------------------------

  pure subroutine rigidBody_update( b, coords )
    class(rigidBody), intent(inout) :: b
    real(rb),         intent(in)    :: coords(3,b%NP)

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
    b%MoI = eigenvalues( inertia )
    b%invMoI = one/b%MoI
    A = transpose(eigenvectors( inertia, b%MoI ))

    ! Compute quaternion:
    b%q = quaternion( A )

    ! Calculate position in the body-fixed frame:
    b%d = matmul( A, b%delta )

  end subroutine rigidBody_update

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

  elemental subroutine rigidBody_rotate( b, dt )
    class(rigidBody), intent(inout) :: b
    real(rb),         intent(in)    :: dt
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
        class(rigidBody), intent(inout) :: b
        integer,          intent(in)    :: k
        real(rb),         intent(in)    :: dt
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
  end subroutine rigidBody_rotate

!---------------------------------------------------------------------------------------------------

  pure function rigidBody_particle_momenta( b ) result( P )
    class(rigidBody), intent(in) :: b
    real(rb)                     :: P(3,b%NP)
    integer  :: k
    real(rb) :: omega(3), At(3,3)
    At = matmul( matrix_Ct(b%q), matrix_B(b%q) )
    omega = matmul(At, half*b%invMoI*matmul( matrix_Bt(b%q), b%pi ))
    forall( k = 1: b%NP)
      P(:,k) = b%M(k)*(b%invMass*b%pcm + cross_product(omega, b%delta(:,k)))
    end forall
  end function rigidBody_particle_momenta

!---------------------------------------------------------------------------------------------------

  function rigidBody_force_torque_virial( b, F ) result( virial )
    class(rigidBody), intent(inout) :: b
    real(rb),         intent(in)    :: F(:,:)
    real(rb)                        :: virial
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
  end function rigidBody_force_torque_virial

!---------------------------------------------------------------------------------------------------

end module ArBee
