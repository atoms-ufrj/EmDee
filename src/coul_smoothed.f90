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

module coul_smoothed_module

use coulModelClass

implicit none

!> Abstract class for model coul_smoothed
!!
!! NOTES: 1) model parameters must be declared individually and tagged with a comment mark "!<>"
!!        2) recognizable parameter types are real(rb) and integer(ib)
!!        3) allocatable one-dimensional arrays (i.e. vectors) are permitted as parameters
!!        4) an integer(ib) scalar parameter - a size - must necessarily succeed every allocatable
!!           parameter or series of equally-sized allocatable parameters.

type, extends(cCoulModel) :: coul_smoothed
  real(rb) :: skin    !<> Smoothing skin width (unit = distance)

  real(rb) :: Rm2 = zero
  real(rb) :: invRm = zero
  real(rb) :: factor = zero

  contains
    procedure :: setup => coul_smoothed_setup
    procedure :: apply_cutoff => coul_smoothed_apply_cutoff
    procedure :: compute => coul_smoothed_compute
    procedure :: energy  => coul_smoothed_energy
    procedure :: virial  => coul_smoothed_virial
end type coul_smoothed

contains

!---------------------------------------------------------------------------------------------------

  subroutine coul_smoothed_setup( model, params, iparams )
    class(coul_smoothed), intent(inout) :: model
    real(rb), optional,   intent(in)    :: params(:)
    integer,  optional,   intent(in)    :: iparams(:)

    ! Model name:
    model%name = "smoothed"

    ! Model parameters:
    model%skin = params(1)

  end subroutine coul_smoothed_setup

!---------------------------------------------------------------------------------------------------
! This subroutine must define pre-computed quantities that depend on the cutoff distance.

  subroutine coul_smoothed_apply_cutoff( model, Rc )
    class(coul_smoothed), intent(inout) :: model
    real(rb),             intent(in)    :: Rc

    model%Rm2 = (Rc - model%skin)**2
    model%invRm = one/(Rc - model%skin)
    model%factor = one/(Rc**2 - model%Rm2)

  end subroutine coul_smoothed_apply_cutoff

!---------------------------------------------------------------------------------------------------
! This subroutine must return the Coulombic energy E(r) and virial W(r) = -r*dE/dr of a pair ij
! whose distance is equal to 1/invR. If the Coulomb model requires a kspace solver, then only the
! real-space, short-range contribution must be computed here.

  subroutine coul_smoothed_compute( model, ECij, WCij, invR, invR2, QiQj )
    class(coul_smoothed), intent(in)  :: model
    real(rb),             intent(out) :: ECij, WCij
    real(rb),             intent(in)  :: invR, invR2, QiQj

    real(rb) :: r2, u, u2, u3, G, WG

    WCij = QiQj*invR
    ECij = WCij

    if (invR < model%invRm) then
      r2 = one/invR2
      u = model%factor*(r2 - model%Rm2)
      u2 = u*u
      u3 = u*u2
      G = 1.0_rb + u3*(15.0_rb*u - 6.0_rb*u2 - 10.0_rb)
      WG = -60.0_rb*u2*(2.0_rb*u - u2 - 1.0_rb)*model%factor*r2
      WCij = WCij*G + ECij*WG
      ECij = ECij*G
    end if

  end subroutine coul_smoothed_compute

!---------------------------------------------------------------------------------------------------
! This subroutine is similar to coulModel_compute, except that only the energy must be computed.

  subroutine coul_smoothed_energy( model, ECij, invR, invR2, QiQj )
    class(coul_smoothed), intent(in)  :: model
    real(rb),             intent(out) :: ECij
    real(rb),             intent(in)  :: invR, invR2, QiQj

    real(rb) :: r2, u, u2, u3, G

    ECij = QiQj*invR

    if (invR < model%invRm) then
      r2 = one/invR2
      u = model%factor*(r2 - model%Rm2)
      u2 = u*u
      u3 = u*u2
      G = 1.0_rb + u3*(15.0_rb*u - 6.0_rb*u2 - 10.0_rb)
      ECij = ECij*G
    end if

  end subroutine coul_smoothed_energy

!---------------------------------------------------------------------------------------------------
! This subroutine is similar to coulModel_compute, except that only the virial must be computed.

  subroutine coul_smoothed_virial( model, WCij, invR, invR2, QiQj )
    class(coul_smoothed), intent(in)  :: model
    real(rb),             intent(out) :: WCij
    real(rb),             intent(in)  :: invR, invR2, QiQj

    real(rb) :: r2, u, u2, u3, G, WG

    WCij = QiQj*invR

    if (invR < model%invRm) then
      r2 = one/invR2
      u = model%factor*(r2 - model%Rm2)
      u2 = u*u
      u3 = u*u2
      G = 1.0_rb + u3*(15.0_rb*u - 6.0_rb*u2 - 10.0_rb)
      WG = -60.0_rb*u2*(2.0_rb*u - u2 - 1.0_rb)*model%factor*r2
      WCij = WCij*(G + WG)
    end if

  end subroutine coul_smoothed_virial

!---------------------------------------------------------------------------------------------------

end module coul_smoothed_module
