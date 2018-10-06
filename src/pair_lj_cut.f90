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

module pair_lj_cut_module

use global
use pairModelClass

implicit none

!> Abstract class for pair model lj_cut
!!
!! NOTES: 1) model parameters must be declared individually and tagged with a comment mark "!<>"
!!        2) recognizable parameter types are real(rb) and integer(ib)
!!        3) allocatable one-dimensional arrays (i.e. vectors) are permitted as parameters
!!        4) an integer(ib) scalar parameter - a size - must necessarily succeed every allocatable
!!           parameter or series of equally-sized allocatable parameters.

type, extends(cPairModel) :: pair_lj_cut
  real(rb) :: epsilon !<> Depth of the potential well
  real(rb) :: sigma   !<> Distance at which the potential is zero

  real(rb) :: eps4, eps24, sigsq
  contains
    procedure :: setup => pair_lj_cut_setup
    procedure :: compute => pair_lj_cut_compute
    procedure :: energy => pair_lj_cut_energy
    procedure :: virial => pair_lj_cut_virial
    procedure :: mix => pair_lj_cut_mix
end type pair_lj_cut

contains

!---------------------------------------------------------------------------------------------------

  subroutine pair_lj_cut_setup( model, params, iparams )
    class(pair_lj_cut), intent(inout) :: model
    real(rb), intent(in), optional :: params(:)
    integer,  intent(in), optional :: iparams(:)

    ! Model name:
    model%name = "lj_cut"

    ! Model parameters:
    model%epsilon = params(1)
    model%sigma = params(2)

    ! Pre-computed quantities:
    model%eps4  = 4.0_rb*model%epsilon
    model%eps24 = 24.0_rb*model%epsilon
    model%sigsq = model%sigma**2

  end subroutine pair_lj_cut_setup

!---------------------------------------------------------------------------------------------------

  subroutine pair_lj_cut_compute( model, Eij, Wij, invR, invR2, QiQj )
    class(pair_lj_cut), intent(in)  :: model
    real(rb),           intent(out) :: Eij, Wij
    real(rb),           intent(in)  :: invR, invR2, QiQj

    real(rb) :: sr2, sr6, sr12

    sr2 = model%sigSq*invR2
    sr6 = sr2*sr2*sr2
    sr12 = sr6*sr6
    Eij = model%eps4*(sr12 - sr6)
    Wij = model%eps24*(sr12 + sr12 - sr6)

  end subroutine pair_lj_cut_compute

!---------------------------------------------------------------------------------------------------

  subroutine pair_lj_cut_energy( model, Eij, invR, invR2, QiQj )
    class(pair_lj_cut), intent(in)  :: model
    real(rb),           intent(out) :: Eij
    real(rb),           intent(in)  :: invR, invR2, QiQj

    real(rb) :: sr2, sr6, sr12

    sr2 = model%sigSq*invR2
    sr6 = sr2*sr2*sr2
    sr12 = sr6*sr6
    Eij = model%eps4*(sr12 - sr6)

  end subroutine pair_lj_cut_energy

!---------------------------------------------------------------------------------------------------

  subroutine pair_lj_cut_virial( model, Wij, invR, invR2, QiQj )
    class(pair_lj_cut), intent(in)  :: model
    real(rb),           intent(out) :: Wij
    real(rb),           intent(in)  :: invR, invR2, QiQj

    real(rb) :: sr2, sr6, sr12

    sr2 = model%sigSq*invR2
    sr6 = sr2*sr2*sr2
    sr12 = sr6*sr6
    Wij = model%eps24*(sr12 + sr12 - sr6)

  end subroutine pair_lj_cut_virial

!---------------------------------------------------------------------------------------------------

  function pair_lj_cut_mix( this, other ) result( mixed )
    class(pair_lj_cut), intent(in) :: this
    class(cPairModel),  intent(in) :: other
    class(cPairModel),  pointer    :: mixed

    select type (other)

      class is (pair_lj_cut)
        allocate(pair_lj_cut :: mixed)
        call mixed % setup( [sqrt(this%epsilon*other%epsilon), half*(this%sigma + other%sigma)] )

      class default
        mixed => null()

    end select

  end function pair_lj_cut_mix

!---------------------------------------------------------------------------------------------------

end module pair_lj_cut_module
