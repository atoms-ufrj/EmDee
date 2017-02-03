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

module pair_lj_sf_module

use global
use pairModelClass

implicit none

!> Abstract class for pair model lj_sf
!!
!! NOTES: 1) model parameters must be declared individually and tagged with a comment mark "!<>"
!!        2) recognizable parameter types are real(rb) and integer(ib)
!!        3) allocatable one-dimensional arrays (i.e. vectors) are permitted as parameters
!!        4) an integer(ib) scalar parameter - a size - must necessarily succeed every allocatable
!!           parameter or series of equally-sized allocatable parameters.

type, extends(cPairModel) :: pair_lj_sf
  real(rb) :: epsilon !<> Depth of the potential well
  real(rb) :: sigma   !<> Distance at which the potential is zero

  real(rb) :: eps4, eps24, sigsq
  contains
    procedure :: setup => pair_lj_sf_setup
    procedure :: compute => pair_lj_sf_compute
    procedure :: virial => pair_lj_sf_virial
    procedure :: mix => pair_lj_sf_mix
end type pair_lj_sf

contains

!---------------------------------------------------------------------------------------------------

  subroutine pair_lj_sf_setup( model, params, iparams )
    class(pair_lj_sf), intent(inout) :: model
    real(rb), intent(in), optional :: params(:)
    integer,  intent(in), optional :: iparams(:)

    ! Model name:
    model%name = "lj_sf"

    ! Model parameters:
    model%epsilon = params(1)
    model%sigma = params(2)

    ! Pre-computed quantities:
    model%eps4 = 4.0_rb*model%epsilon
    model%eps24 = 24.0_rb*model%epsilon
    model%sigsq = model%sigma**2

    ! Activate shifted-force status:
    model%shifted_force = .true.

    ! Attest that invR is evaluated during model computations:
    model%noInvR = .false.
    model%noInvR_virial = .false.

  end subroutine pair_lj_sf_setup

!---------------------------------------------------------------------------------------------------

  subroutine pair_lj_sf_compute( model, Eij, Wij, invR, invR2 )
    class(pair_lj_sf), intent(in)  :: model
    real(rb),          intent(out) :: Eij, Wij, invR
    real(rb),          intent(in)  :: invR2

    real(rb) :: sr2, sr6, sr12, rFc

    sr2 = model%sigSq*invR2
    sr6 = sr2*sr2*sr2
    sr12 = sr6*sr6
    invR = sqrt(invR2)
    rFc = model%fshift/invR
    Eij = model%eps4*(sr12 - sr6) + model%eshift + rFc
    Wij = model%eps24*(sr12 + sr12 - sr6) - rFc

  end subroutine pair_lj_sf_compute

!---------------------------------------------------------------------------------------------------

  subroutine pair_lj_sf_virial( model, Wij, invR, invR2 )
    class(pair_lj_sf), intent(in)  :: model
    real(rb),          intent(out) :: Wij, invR
    real(rb),          intent(in)  :: invR2

    real(rb) :: sr2, sr6, sr12

    sr2 = model%sigSq*invR2
    sr6 = sr2*sr2*sr2
    sr12 = sr6*sr6
    invR = sqrt(invR2)
    Wij = model%eps24*(sr12 + sr12 - sr6) - model%fshift/invR

  end subroutine pair_lj_sf_virial

!---------------------------------------------------------------------------------------------------

  function pair_lj_sf_mix( this, other ) result( mixed )
    class(pair_lj_sf), intent(in) :: this
    class(cPairModel), intent(in) :: other
    class(cPairModel), pointer    :: mixed

    select type (other)

      class is (pair_lj_sf)
        allocate(pair_lj_sf :: mixed)
        call mixed % setup( [sqrt(this%epsilon*other%epsilon), half*(this%sigma + other%sigma)] )

      class default
        mixed => null()

    end select

  end function pair_lj_sf_mix

!---------------------------------------------------------------------------------------------------

end module pair_lj_sf_module
