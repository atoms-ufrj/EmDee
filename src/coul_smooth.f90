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

module coul_smooth_module

use coulModelClass
use math, only : uerfc

implicit none

!> Abstract class for model coul_smooth
!!
!! NOTES: 1) model parameters must be declared individually and tagged with a comment mark "!<>"
!!        2) recognizable parameter types are real(rb) and integer(ib)
!!        3) allocatable one-dimensional arrays (i.e. vectors) are permitted as parameters
!!        4) an integer(ib) scalar parameter - a size - must necessarily succeed every allocatable
!!           parameter or series of equally-sized allocatable parameters.

type, extends(cCoulModel) :: coul_smooth
  real(rb) :: Damp    !<> Damping parameter (unit = 1/Distance)
  real(rb) :: Rm      !<> Distance at which smoothing starts
  real(rb) :: Rc      !<> Cut-off distance

  real(rb) :: beta, Rm2, invRm, factor

  contains
    procedure :: setup => coul_smooth_setup
    procedure :: kspace_setup => coul_smooth_kspace_setup
    procedure :: compute => coul_smooth_compute
    procedure :: energy  => coul_smooth_energy
    procedure :: virial  => coul_smooth_virial
end type coul_smooth

contains

!---------------------------------------------------------------------------------------------------

  subroutine coul_smooth_setup( model, params, iparams )
    class(coul_smooth), intent(inout) :: model
    real(rb), optional, intent(in)    :: params(:)
    integer,  optional, intent(in)    :: iparams(:)

    ! Model name:
    model%name = "smooth"

    ! Model parameters:
    model%Damp = params(1)
    model%Rm = params(2)
    model%Rc = params(3)

    ! Pre-computed quantities:
    model%Rm2 = model%Rm**2
    model%invRm = one/model%Rm
    model%factor = one/(model%Rc**2 - model%Rm2)
    model%alpha = model%Damp
    model%beta = two*model%alpha/sqrt(Pi)

  end subroutine coul_smooth_setup

!---------------------------------------------------------------------------------------------------
! For a Coulomb model that requires a kspace solver, this subroutine must be used to store the
! values of parameters that depend on the Ewald splitting constant 'alpha'.

  subroutine coul_smooth_kspace_setup( model, alpha )
    class(coul_smooth), intent(inout) :: model
    real(rb),           intent(in)    :: alpha

    model%alpha = model%Damp
    model%beta = two*alpha/sqrt(Pi)

  end subroutine coul_smooth_kspace_setup

!---------------------------------------------------------------------------------------------------
! This subroutine must return the Coulombic energy E(r) and virial W(r) = -r*dE/dr of a pair ij
! whose distance is equal to 1/invR. If the Coulomb model requires a kspace solver, then only the
! real-space, short-range contribution must be computed here.

  subroutine coul_smooth_compute( model, ECij, WCij, invR, invR2, QiQj )
    class(coul_smooth), intent(in)  :: model
    real(rb),           intent(out) :: ECij, WCij
    real(rb),           intent(in)  :: invR, invR2, QiQj

    real(rb) :: x, expmx2, r2, u, u2, u3, G, WG

    x = model%alpha/invR
    expmx2 = exp(-x*x)
    ECij = QiQj*uerfc(x,expmx2)*invR
    WCij = ECij + QiQj*model%beta*expmx2

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

  end subroutine coul_smooth_compute

!---------------------------------------------------------------------------------------------------
! This subroutine is similar to coulModel_compute, except that only the energy must be computed.

  subroutine coul_smooth_energy( model, ECij, invR, invR2, QiQj )
    class(coul_smooth), intent(in)  :: model
    real(rb),           intent(out) :: ECij
    real(rb),           intent(in)  :: invR, invR2, QiQj

    real(rb) :: x, r2, u, u2, u3, G

    x = model%alpha/invR
    ECij = QiQj*uerfc(x,exp(-x*x))*invR

    if (invR < model%invRm) then
      r2 = one/invR2
      u = model%factor*(r2 - model%Rm2)
      u2 = u*u
      u3 = u*u2
      G = 1.0_rb + u3*(15.0_rb*u - 6.0_rb*u2 - 10.0_rb)
      ECij = ECij*G
    end if

  end subroutine coul_smooth_energy

!---------------------------------------------------------------------------------------------------
! This subroutine is similar to coulModel_compute, except that only the virial must be computed.

  subroutine coul_smooth_virial( model, WCij, invR, invR2, QiQj )
    class(coul_smooth), intent(in)  :: model
    real(rb),           intent(out) :: WCij
    real(rb),           intent(in)  :: invR, invR2, QiQj

    real(rb) :: x, ECij, expmx2, r2, u, u2, u3, G, WG

    x = model%alpha/invR
    expmx2 = exp(-x*x)
    ECij = QiQj*uerfc(x,expmx2)*invR
    WCij = ECij + QiQj*model%beta*expmx2

    if (invR < model%invRm) then
      r2 = one/invR2
      u = model%factor*(r2 - model%Rm2)
      u2 = u*u
      u3 = u*u2
      G = 1.0_rb + u3*(15.0_rb*u - 6.0_rb*u2 - 10.0_rb)
      WG = -60.0_rb*u2*(2.0_rb*u - u2 - 1.0_rb)*model%factor*r2
      WCij = WCij*G + ECij*WG
    end if

  end subroutine coul_smooth_virial

!---------------------------------------------------------------------------------------------------

end module coul_smooth_module
