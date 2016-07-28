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
  real(rb) :: pcm(3)    ! Center-of-mass momentum
  real(rb) :: q(4)      ! Unit quaternion of orientation
  real(rb) :: pi(4)     ! Quaternion momentum

  integer,  allocatable :: index(:)
  real(rb), allocatable :: M(:)
  real(rb), allocatable :: d(:,:)

  real(rb) :: invMass    ! Inverse of body mass
  real(rb) :: invMoI(3)  ! Inverses of principal moments of inertia

  contains

    procedure :: setup => rigidBody_setup
    procedure :: update => rigidBody_update

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

  subroutine rigidBody_setup( body, coords, L, masses, indexes )
    class(rigidBody), intent(inout) :: body
    real(rb),         intent(inout) :: coords(:,:)
    real(rb),         intent(in)    :: L(3), masses(:)
    integer(ib),      intent(in)    :: indexes(:)

    integer  :: k, i
    real(rb) :: r1(3)

    body%NP = size(indexes)
    allocate( body%index(body%NP), body%M(body%NP), body%d(3,body%NP) )
    body%index = indexes
    body%M = masses(indexes)
    body%mass = sum(body%M)
    body%invMass = one/body%mass
    r1 = coords(:,indexes(1))
    do k = 2, body%NP
      i = indexes(k)
      coords(:,i) = coords(:,i) - L*anint((coords(:,i) - r1)/L)
    end do
    body%pcm = zero
    body%pi  = zero
    call body % update( coords )
  end subroutine rigidBody_setup

!---------------------------------------------------------------------------------------------------

  subroutine rigidBody_update( body, coords )
    class(rigidBody), intent(inout) :: body
    real(rb),         intent(in)    :: coords(:,:)

    integer  :: x
    real(rb) :: delta(3,body%NP), inertia(3,3), At(3,3)

    ! Compute center-of-mass position and space-frame particle coordinates:
    delta = coords(:,body%index)
    forall (x=1:3) body%rcm(x) = sum(body%M*delta(x,:))*body%invMass
    forall (x=1:3) delta(x,:) = delta(x,:) - body%rcm(x)

    ! Compute upper-triangular inertia tensor:
    inertia(1,1) = sum(body%M*(delta(2,:)**2 + delta(3,:)**2))
    inertia(2,2) = sum(body%M*(delta(1,:)**2 + delta(3,:)**2))
    inertia(3,3) = sum(body%M*(delta(1,:)**2 + delta(2,:)**2))
    inertia(1,2) = -sum(body%M*delta(1,:)*delta(2,:))
    inertia(1,3) = -sum(body%M*delta(1,:)*delta(3,:))
    inertia(2,3) = -sum(body%M*delta(2,:)*delta(3,:))

    ! Diagonalize the inertia tensor:
    body%MoI = eigenvalues( inertia )
    body%invMoI = one/body%MoI
    At = eigenvectors( inertia, body%MoI )

    ! Compute quaternion:
    body%q = quaternion( transpose(At) )

    ! Calculate position in the body-fixed frame:
    body%d = matmul( delta, At )

  end subroutine rigidBody_update

!---------------------------------------------------------------------------------------------------

end module ArBee
