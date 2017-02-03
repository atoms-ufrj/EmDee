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

module angle_harmonic_module

use global
use angleModelClass

implicit none

!> Abstract class for angle model harmonic
!!
!! NOTES: 1) model parameters must be declared individually and tagged with a comment mark "!<>"
!!        2) recognizable parameter types are real(rb) and integer(ib)
!!        3) allocatable one-dimensional arrays (i.e. vectors) are permitted as parameters
!!        4) an integer(ib) scalar parameter - a size - must necessarily succeed every allocatable
!!           parameter or series of equally-sized allocatable parameters.

type, extends(cAngleModel) :: angle_harmonic
  real(rb) :: k       !<> Force constant
  real(rb) :: theta0  !<> Equilibrium angle

  real(rb) :: minus_k, half_k
  contains
    procedure :: setup => angle_harmonic_setup
    procedure :: compute => angle_harmonic_compute
end type angle_harmonic

contains

!---------------------------------------------------------------------------------------------------

  subroutine angle_harmonic_setup( model, params, iparams )
    class(angle_harmonic), intent(inout) :: model
    real(rb), intent(in), optional :: params(:)
    integer,  intent(in), optional :: iparams(:)

    ! Model name:
    model%name = "harmonic"

    ! Model parameters:
    model%k = params(1)
    model%theta0 = params(2)

    ! Pre-computed quantities:
    model%minus_k = -model%k
    model%half_k = half*model%k

  end subroutine angle_harmonic_setup

!---------------------------------------------------------------------------------------------------

  subroutine angle_harmonic_compute( model, Ea, Fa, theta )
    class(angle_harmonic), intent(in)  :: model
    real(rb),              intent(out) :: Ea, Fa
    real(rb),              intent(in)  :: theta

    real(rb) :: delta
    delta = theta - model%theta0
    Ea = model%half_k*delta**2
    Fa = model%minus_k*delta

  end subroutine angle_harmonic_compute

!---------------------------------------------------------------------------------------------------

end module angle_harmonic_module
