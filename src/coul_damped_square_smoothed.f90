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

module coul_damped_square_smoothed_module

use coulModelClass
use math, only : uerfc

implicit none

!> Abstract class for model coul_damped_square_smoothed
!!
!! NOTES: 1) model parameters must be declared individually and tagged with a comment mark "!<>"
!!        2) recognizable parameter types are real(rb) and integer(ib)
!!        3) allocatable one-dimensional arrays (i.e. vectors) are permitted as parameters
!!        4) an integer(ib) scalar parameter - a size - must necessarily succeed every allocatable
!!           parameter or series of equally-sized allocatable parameters.

type, extends(cCoulModel) :: coul_damped_square_smoothed
  real(rb) :: damp    !<> Damping parameter (unit = 1/distance)
  real(rb) :: skinWidth    !<> Smoothing skinWidth width (unit = distance)

  real(rb) :: beta = zero
  real(rb) :: Rm2 = zero
  real(rb) :: invRm = zero

  contains
    procedure :: setup => coul_damped_square_smoothed_setup
    procedure :: apply_cutoff => coul_damped_square_smoothed_apply_cutoff
    procedure :: compute => coul_damped_square_smoothed_compute
    procedure :: energy  => coul_damped_square_smoothed_energy
    procedure :: virial  => coul_damped_square_smoothed_virial
end type coul_damped_square_smoothed

contains

!---------------------------------------------------------------------------------------------------

  subroutine coul_damped_square_smoothed_setup( model, params, iparams )
    class(coul_damped_square_smoothed), intent(inout) :: model
    real(rb), optional,          intent(in)    :: params(:)
    integer,  optional,          intent(in)    :: iparams(:)

    ! Model name:
    model%name = "damped_smoothed"

    ! Model parameters:
    model%Damp = params(1)
    model%skinWidth = params(2)

    ! Pre-computed quantities:
    model%alpha = model%damp
    model%beta = two*model%alpha/sqrt(Pi)

  end subroutine coul_damped_square_smoothed_setup

!---------------------------------------------------------------------------------------------------
! This subroutine must define pre-computed quantities that depend on the cutoff distance.

  subroutine coul_damped_square_smoothed_apply_cutoff( model, Rc )
    class(coul_damped_square_smoothed), intent(inout) :: model
    real(rb),                    intent(in)    :: Rc

    model%Rm2 = (Rc - model%skinWidth)**2
    model%invRm = one/(Rc - model%skinWidth)
    model%factor = one/(Rc**2 - model%Rm2)

  end subroutine coul_damped_square_smoothed_apply_cutoff

!---------------------------------------------------------------------------------------------------
! This subroutine must return the Coulombic energy E(r) and virial W(r) = -r*dE/dr of a pair ij
! whose distance is equal to 1/invR. If the Coulomb model requires a kspace solver, then only the
! real-space, short-range contribution must be computed here.

  subroutine coul_damped_square_smoothed_compute( model, Eij, Wij, invR, invR2 )
    class(coul_damped_square_smoothed), intent(in)  :: model
    real(rb),                    intent(out) :: Eij, Wij
    real(rb),                    intent(in)  :: invR, invR2

    real(rb) :: x, expmx2, r2, u, u2, u3, G, WG

    x = model%alpha/invR
    expmx2 = exp(-x*x)
    Eij = uerfc(x,expmx2)*invR
    Wij = Eij + model%beta*expmx2

    if (invR < model%invRm) then
      r2 = one/invR2
      u = model%factor*(r2 - model%Rm2)
      u2 = u*u
      u3 = u*u2
      G = 1.0_rb + u3*(15.0_rb*u - 6.0_rb*u2 - 10.0_rb)
      WG = -60.0_rb*u2*(2.0_rb*u - u2 - 1.0_rb)*model%factor*r2
      Wij = Wij*G + Eij*WG
      Eij = Eij*G
    end if

  end subroutine coul_damped_square_smoothed_compute

!---------------------------------------------------------------------------------------------------
! This subroutine is similar to coulModel_compute, except that only the energy must be computed.

  subroutine coul_damped_square_smoothed_energy( model, Eij, invR, invR2 )
    class(coul_damped_square_smoothed), intent(in)  :: model
    real(rb),                    intent(out) :: Eij
    real(rb),                    intent(in)  :: invR, invR2

    real(rb) :: x, r2, u, u2, u3, G

    x = model%alpha/invR
    Eij = uerfc(x,exp(-x*x))*invR

    if (invR < model%invRm) then
      r2 = one/invR2
      u = model%factor*(r2 - model%Rm2)
      u2 = u*u
      u3 = u*u2
      G = 1.0_rb + u3*(15.0_rb*u - 6.0_rb*u2 - 10.0_rb)
      Eij = Eij*G
    end if

  end subroutine coul_damped_square_smoothed_energy

!---------------------------------------------------------------------------------------------------
! This subroutine is similar to coulModel_compute, except that only the virial must be computed.

  subroutine coul_damped_square_smoothed_virial( model, Wij, invR, invR2 )
    class(coul_damped_square_smoothed), intent(in)  :: model
    real(rb),                    intent(out) :: Wij
    real(rb),                    intent(in)  :: invR, invR2

    real(rb) :: x, Eij, expmx2, r2, u, u2, u3, G, WG

    x = model%alpha/invR
    expmx2 = exp(-x*x)
    Eij = uerfc(x,expmx2)*invR
    Wij = Eij + model%beta*expmx2

    if (invR < model%invRm) then
      r2 = one/invR2
      u = model%factor*(r2 - model%Rm2)
      u2 = u*u
      u3 = u*u2
      G = 1.0_rb + u3*(15.0_rb*u - 6.0_rb*u2 - 10.0_rb)
      WG = -60.0_rb*u2*(2.0_rb*u - u2 - 1.0_rb)*model%factor*r2
      Wij = Wij*G + Eij*WG
    end if

  end subroutine coul_damped_square_smoothed_virial

!---------------------------------------------------------------------------------------------------

end module coul_damped_square_smoothed_module
