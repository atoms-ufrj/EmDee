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

module coul_long_module

use coulModelClass
use math, only : uerfc

implicit none

!> Abstract class for model coul_long
!!
!! NOTES: 1) model parameters must be declared individually and tagged with a comment mark "!<>"
!!        2) recognizable parameter types are real(rb) and integer(ib)
!!        3) allocatable one-dimensional arrays (i.e. vectors) are permitted as parameters
!!        4) an integer(ib) scalar parameter - a size - must necessarily succeed every allocatable
!!           parameter or series of equally-sized allocatable parameters.

type, extends(cCoulModel) :: coul_long
  real(rb) :: beta
  contains
    procedure :: setup => coul_long_setup
    procedure :: kspace_setup => coul_long_kspace_setup
    procedure :: compute => coul_long_compute
    procedure :: energy  => coul_long_energy
    procedure :: virial  => coul_long_virial
end type coul_long

contains

!---------------------------------------------------------------------------------------------------

  subroutine coul_long_setup( model, params, iparams )
    class(coul_long),   intent(inout) :: model
    real(rb), optional, intent(in)    :: params(:)
    integer,  optional, intent(in)    :: iparams(:)

    ! Model name:
    model%name = "long"

    ! Attest that the model requires kspace:
    model%requires_kspace = .true.

  end subroutine coul_long_setup

!---------------------------------------------------------------------------------------------------
! For a Coulomb model that requires a kspace solver, this subroutine must be used to store the
! values of parameters that depend on the Ewald splitting constant 'alpha'.

  subroutine coul_long_kspace_setup( model, alpha )
    class(coul_long), intent(inout) :: model
    real(rb),         intent(in)    :: alpha

    model%alpha = alpha
    model%beta = two*alpha/sqrt(Pi)

  end subroutine coul_long_kspace_setup

!---------------------------------------------------------------------------------------------------
! This subroutine must return the Coulombic energy E(r) and virial W(r) = -r*dE/dr of a pair ij
! whose distance is equal to 1/invR. If the Coulomb model requires a kspace solver, then only the
! real-space, short-range contribution must be computed here.

  subroutine coul_long_compute( model, Eij, Wij, invR, invR2 )
    class(coul_long), intent(in)  :: model
    real(rb),         intent(out) :: Eij, Wij
    real(rb),         intent(in)  :: invR, invR2

    real(rb) :: x, expmx2

    x = model%alpha/invR
    expmx2 = exp(-x*x)
    Eij = uerfc(x,expmx2)*invR
    Wij = Eij + model%beta*expmx2

  end subroutine coul_long_compute

!---------------------------------------------------------------------------------------------------
! This subroutine is similar to coulModel_compute, except that only the energy must be computed.

  subroutine coul_long_energy( model, Eij, invR, invR2 )
    class(coul_long), intent(in)  :: model
    real(rb),         intent(out) :: Eij
    real(rb),         intent(in)  :: invR, invR2

    real(rb) :: x

    x = model%alpha/invR
    Eij = uerfc(x,exp(-x*x))*invR

  end subroutine coul_long_energy

!---------------------------------------------------------------------------------------------------
! This subroutine is similar to coulModel_compute, except that only the virial must be computed.

  subroutine coul_long_virial( model, Wij, invR, invR2 )
    class(coul_long), intent(in)  :: model
    real(rb),         intent(out) :: Wij
    real(rb),         intent(in)  :: invR, invR2

    real(rb) :: x, expmx2

    x = model%alpha/invR
    expmx2 = exp(-x*x)
    Wij = uerfc(x,expmx2)*invR + model%beta*expmx2

  end subroutine coul_long_virial

!---------------------------------------------------------------------------------------------------

end module coul_long_module
