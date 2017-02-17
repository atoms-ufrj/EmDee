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

module pair_mie_cut_module

use global
use pairModelClass

implicit none

!> Abstract class for pair model mie
!!
!! NOTES: 1) model parameters must be declared individually and tagged with a comment mark "!<>"
!!        2) recognizable parameter types are real(rb) and integer(ib)
!!        3) allocatable one-dimensional arrays (i.e. vectors) are permitted as parameters
!!        4) an integer(ib) scalar parameter - a size - must necessarily succeed every allocatable
!!           parameter or series of equally-sized allocatable parameters.

type, extends(cPairModel) :: pair_mie_cut
  real(rb) :: epsilon !<> Depth of the potential well
  real(rb) :: sigma   !<> Distance at which the potential is zero
  real(rb) :: n !<> Softness/hardness of the repulsion
  real(rb) :: m !<>  Range of the attraction of the intermolecular potential
  real(rb) :: miem,mien
  contains
    procedure :: setup => pair_mie_cut_setup
    procedure :: compute => pair_mie_cut_compute
    procedure :: virial => pair_mie_cut_virial
    procedure :: mix => pair_mie_cut_mix
end type pair_mie_cut

contains

!---------------------------------------------------------------------------------------------------

  subroutine pair_mie_cut_setup( model, params, iparams )
    class(pair_mie_cut), intent(inout) :: model
    real(rb), intent(in), optional :: params(:)
    integer,  intent(in), optional :: iparams(:)
    real(rb) :: coeff
    ! Model name:
    model%name = "mie_cut"

    ! Model parameters:
    model%epsilon = params(1)
    model%sigma = params(2)
    model%n = params(3)
    model%m = params (4)
    ! Pre-computed quantities:
    coeff  = (model%n / (model%n - model%m)) * (model%n / model%m) ** (model%m / (model%n - model%m)) * model%epsilon
    model%mien = coeff * model%sigma ** model%n
    model%miem = coeff * model%sigma ** model%m

  end subroutine pair_mie_cut_setup

!---------------------------------------------------------------------------------------------------

  subroutine pair_mie_cut_compute( model, Eij, Wij, invR, invR2 )
    class(pair_mie_cut), intent(in)  :: model
    real(rb),           intent(out) :: Eij, Wij, invR
    real(rb),           intent(in)  :: invR2

    real(rb) :: rn, rm

    rm = invr2 ** (model%m / 2.0_rb)
    rn = rm ** (model%n / model%m)

    Eij = model%mien * rn - model%miem * rm 
    Wij = model%n * model%mien * rn - model%m * model%miem * rm

  end subroutine pair_mie_cut_compute

!---------------------------------------------------------------------------------------------------

  subroutine pair_mie_cut_virial( model, Wij, invR, invR2 )
    class(pair_mie_cut), intent(in)  :: model
    real(rb),           intent(out) :: Wij, invR
    real(rb),           intent(in)  :: invR2

    real(rb) :: rn, rm

    rm = invr2 ** (model%m / 2.0_rb)
    rn = rm ** (model%n / model%m)

    Wij = model%n * model%mien * rn - model%m * model%miem * rm

  end subroutine pair_mie_cut_virial

!---------------------------------------------------------------------------------------------------

  function pair_mie_cut_mix( this, other ) result( mixed )
    class(pair_mie_cut), intent(in) :: this
    class(cPairModel),  intent(in) :: other
    class(cPairModel),  pointer    :: mixed

    select type (other)

      class is (pair_mie_cut)
        allocate(pair_mie_cut :: mixed)
        call mixed % setup( [sqrt(this%epsilon*other%epsilon)*  &
                            sqrt(this%sigma**3.0_rb*other%sigma**3.0_rb) & 
                            /(half*(this%sigma + other%sigma)), &
                             half*(this%sigma + other%sigma),  &
                             3.0_rb + sqrt((this%n - 3.0_rb)*(other%n - 3.0_rb)),  &
                             3.0_rb + sqrt((this%m - 3.0_rb)*(other%m - 3.0_rb))] )
      class default
        mixed => null()

    end select

  end function pair_mie_cut_mix

!---------------------------------------------------------------------------------------------------

end module pair_mie_cut_module
