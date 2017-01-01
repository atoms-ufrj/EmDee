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

module pair_lj_cut_coul_sf_module

use pairModelClass
use pair_lj_cut_module

!> Abstract class for pair model lj_cut_coul_sf
!! NOTES: 1) model parameters must be declared individually and tagged with a comment mark "!<>"
!!        2) allocatable parameters are permitted only for rank=1
!!        3) a series of rank-1 allocatable parameters must be succeeded by an integer parameter,
!!           which will contain their (common) size after allocation
type, extends(cPairModel) :: pair_lj_cut_coul_sf
  real(rb) :: epsilon !<> Lennard-Jones parameter epsilon
  real(rb) :: sigma   !<> Lennard-Jones parameter sigma

  real(rb) :: eps4, eps24, sigsq
  contains
    procedure :: setup => pair_lj_cut_coul_sf_setup
    procedure :: compute => pair_lj_cut_coul_sf_compute
    procedure :: virial => pair_lj_cut_coul_sf_virial
    procedure :: mix => pair_lj_cut_coul_sf_mix
end type pair_lj_cut_coul_sf

contains

!---------------------------------------------------------------------------------------------------

  subroutine pair_lj_cut_coul_sf_setup( model, params, iparams )
    class(pair_lj_cut_coul_sf), intent(inout) :: model
    real(rb), intent(in), optional :: params(:)
    integer,  intent(in), optional :: iparams(:)

    ! Model name:
    model%name = "lj_cut_coul_sf"

    ! Model parameters:
    model%epsilon = params(1)
    model%sigma = params(2)

    ! Pre-computed quantities:
    model%eps4 = 4.0_rb*model%epsilon
    model%eps24 = 24.0_rb*model%epsilon
    model%sigsq = model%sigma**2

    ! Mark active contributions:
    model%vdw = .true.
    model%coulomb = .true.

    ! Activate shifted-force status:
    model%shifted_force_coul = .true.

  end subroutine pair_lj_cut_coul_sf_setup

!---------------------------------------------------------------------------------------------------

  subroutine pair_lj_cut_coul_sf_compute( model, Eij, Wij, invR2, Qi, Qj )
    class(pair_lj_cut_coul_sf), intent(in)  :: model
    real(rb),                   intent(out) :: Eij, Wij
    real(rb),                   intent(in)  :: invR2, Qi, Qj

    real(rb) :: sr2, sr6, sr12, rFc, invR, QiQj, QiQjbyR

    sr2 = model%sigSq*invR2
    sr6 = sr2*sr2*sr2
    sr12 = sr6*sr6
    invR = sqrt(invR2)
    QiQj = Qi*Qj
    QiQjbyR = QiQj*invR
    rFc = QiQj*model%fshift_coul/invR
    Eij = model%eps4*(sr12 - sr6) + QiQjbyR + QiQj*model%eshift_coul + rFc
    Wij = model%eps24*(sr12 + sr12 - sr6) + QiQjbyR - rFc

  end subroutine pair_lj_cut_coul_sf_compute

!---------------------------------------------------------------------------------------------------

  function pair_lj_cut_coul_sf_virial( model, invR2, Qi, Qj ) result( Wij )
    class(pair_lj_cut_coul_sf), intent(in) :: model
    real(rb),                   intent(in) :: invR2, Qi, Qj
    real(rb)                               :: Wij

    real(rb) :: sr2, sr6, sr12, invR

    sr2 = model%sigSq*invR2
    sr6 = sr2*sr2*sr2
    sr12 = sr6*sr6
    invR = sqrt(invR2)
    Wij = model%eps24*(sr12 + sr12 - sr6) + Qi*Qj*(invR - model%fshift_coul/invR)

  end function pair_lj_cut_coul_sf_virial

!---------------------------------------------------------------------------------------------------

  function pair_lj_cut_coul_sf_mix( this, other ) result( mixed )
    class(pair_lj_cut_coul_sf),    intent(in) :: this
    class(cPairModel), intent(in) :: other
    class(cPairModel), pointer :: mixed

    select type (other)
      class is (pair_lj_cut_coul_sf)
        allocate(pair_lj_cut_coul_sf :: mixed)
        call mixed % setup( [sqrt(this%epsilon*other%epsilon), half*(this%sigma + other%sigma)] )
      class is (pair_lj_cut)
        allocate(pair_lj_cut :: mixed)
        call mixed % setup( [sqrt(this%epsilon*other%epsilon), half*(this%sigma + other%sigma)] )
      class default
        mixed => null()
    end select

  end function pair_lj_cut_coul_sf_mix

!---------------------------------------------------------------------------------------------------

end module pair_lj_cut_coul_sf_module
