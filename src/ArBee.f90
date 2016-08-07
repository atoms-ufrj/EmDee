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
  real(rb) :: mass      ! Total body mass
  real(rb) :: MoI(3)    ! Principal moments of inertia
  real(rb) :: rcm(3)    ! Center-of-mass position
  real(rb) :: pcm(3)    ! Center-of-mass momentum vector
  real(rb) :: q(0:3)    ! Unit quaternion of orientation
  real(rb) :: pi(4)     ! Quaternion momentum

  integer,  allocatable :: index(:)
  real(rb), allocatable :: M(:)
  real(rb), allocatable :: d(:,:)
  real(rb), allocatable :: delta(:,:)
  real(rb), allocatable :: P(:,:)

  real(rb) :: invMass    ! Inverse of body mass
  real(rb) :: invMoI(3)  ! Inverses of principal moments of inertia

  contains

    procedure :: setup => rigidBody_setup
    procedure :: update => rigidBody_update
    procedure :: B => rigidBody_B
    procedure :: set_momenta => rigidBody_set_momenta

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
    allocate( b%index(b%NP), b%M(b%NP), b%d(3,b%NP), b%delta(3,b%NP), b%P(3,b%NP) )
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
    b%P = zero

  end subroutine rigidBody_setup

!---------------------------------------------------------------------------------------------------

  pure subroutine rigidBody_update( b, coords )
    class(rigidBody), intent(inout) :: b
    real(rb),         intent(in)    :: coords(3,b%NP)

    integer  :: x
    real(rb) :: inertia(3,3), At(3,3)

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
    At = eigenvectors( inertia, b%MoI )

    ! Compute quaternion:
    b%q = quaternion( transpose(At) )

    ! Calculate position in the body-fixed frame:
    b%d = matmul( b%delta, At )

  end subroutine rigidBody_update

!---------------------------------------------------------------------------------------------------

  pure function rigidBody_B( body ) result( B )
    class(rigidBody), intent(in) :: body
    real(rb)                     :: B(4,3)
    associate (q => body%q)
      B = reshape( [-q(1),  q(0),  q(3), -q(2),  &
                    -q(2), -q(3),  q(0),  q(1),  &
                    -q(3),  q(2), -q(1),  q(0)], [4,3] )
    end associate
  end function rigidBody_B

!---------------------------------------------------------------------------------------------------

  subroutine rigidBody_set_momenta( b, pcm, omega )
    class(rigidBody), intent(inout) :: b
    real(rb),         intent(in)    :: pcm(3), omega(3)
    integer :: k
    b%pcm = pcm
    b%pi = matmul(b%B(), two*b%MoI*omega)
    do k = 1, b%NP
      b%P(:,k) = b%M(k)*(b%invMass*b%pcm + cross_product(omega, b%delta(:,k)))
    end do
  end subroutine rigidBody_set_momenta

!---------------------------------------------------------------------------------------------------

end module ArBee
