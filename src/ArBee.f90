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

type, private :: tParticle
  real(rb) :: index   ! System-based index
  real(rb) :: mass    ! Mass
  real(rb) :: d(3)    ! Body-fixed internal coordinates
end type tParticle

type rigidBody
  integer  :: NP      ! Number of particles
  real(rb) :: mass    ! Total mass
  real(rb) :: MoI(3)  ! Principal moments of inertia
  real(rb) :: rcm(3)  ! Center-of-mass position
  real(rb) :: pcm(3)  ! Center-of-mass momentum
  real(rb) :: q(4)    ! Unit quaternion of orientation
  real(rb) :: pi(4)   ! Quaternion momentum
  type(tParticle), allocatable :: particle(:)

  real(rb) :: invMass
  real(rb) :: invMoI(3)

  contains
    procedure :: setup => rigidBody_setup
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

  subroutine rigidBody_setup( me, indexes, coords, masses )
    class(rigidBody), intent(inout) :: me
    integer(ib),      intent(in)    :: indexes(:)
    real(rb),         intent(in)    :: coords(3,size(indexes)), masses(size(indexes))

    integer  :: dir, i
    real(rb) :: delta(3,size(indexes)), inertia(3,3), A(3,3)

    ! Allocate particles:
    me%NP = size(indexes)
    allocate( me%particle(me%NP) )

    ! Save particle masses and positions:
    me%particle%mass = masses

    ! Compute total mass and center-of-mass position:
    me%mass = sum(masses)
    me%invMass = one/me%mass
    forall (dir=1:3) me%rcm(dir) = sum(masses*coords(dir,:))*me%invMass

    ! Compute inertia tensor:
    forall (dir=1:3) delta(dir,:) = coords(dir,:) - me%rcm(dir)
    inertia = zero
    do i = 1, me%NP
      call add_inertia( inertia, me%particle(i)%mass, delta(:,i) )
    end do
    inertia(2,1) = inertia(1,2)
    inertia(3,1) = inertia(1,3)
    inertia(3,2) = inertia(2,3)

    ! Diagonalize the inertia tensor:
    me%MoI = eigenvalues( inertia )
    me%invMoI = one/me%MoI
    A = transpose(eigenvectors( inertia, me%MoI ))

    ! Compute quaternion:
    me%q = quaternion( A )

    ! Calculate position in the body-fixed frame:
    forall (i=1:me%NP) me%particle(i)%d = matmul( A, delta(:,i) )

    ! Zero center-of-mass and quaternion-conjugated momenta:
    me%pcm = zero
    me%pi  = zero

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine add_inertia( inertia, mass, delta )
        real(rb), intent(inout) :: inertia(3,3)
        real(rb), intent(in)    :: mass, delta(3)
        inertia(1,1) = inertia(1,1) + mass*(delta(2)**2 + delta(3)**2)
        inertia(2,2) = inertia(2,2) + mass*(delta(1)**2 + delta(3)**2)
        inertia(3,3) = inertia(3,3) + mass*(delta(1)**2 + delta(2)**2)
        inertia(1,2) = inertia(1,2) - mass*delta(1)*delta(2)
        inertia(1,3) = inertia(1,3) - mass*delta(1)*delta(3)
        inertia(2,3) = inertia(2,3) - mass*delta(2)*delta(3)
      end subroutine add_inertia
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine rigidBody_setup

!---------------------------------------------------------------------------------------------------

end module ArBee
