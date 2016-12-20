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

module pair_softcore_sf_coul_sf_module

use global
use pairModelClass
!use pair_softcore_sf_module
use pair_lj_sf_module
use pair_lj_sf_coul_sf_module

real(rb), parameter, private :: alpha      = 0.5_rb, &
                                exponent_n = 1.0_rb, &
                                exponent_p = 1.0_rb

!> Abstract class for pair model softcore_sf_coul_sf
!! NOTES: 1) model parameters must be declared individually and tagged with a comment mark "!<>"
!!        2) allocatable parameters are permitted only for rank=1
!!        3) a series of rank-1 allocatable parameters must be succeeded by an integer parameter,
!!           which will contain their (common) size after allocation
type, extends(cPairModel) :: pair_softcore_sf_coul_sf
  real(rb) :: epsilon !<> Lennard-Jones parameter epsilon
  real(rb) :: sigma   !<> Lennard-Jones parameter sigma
  real(rb) :: lambda  !<> Coupling parameter

  real(rb) :: prefactor, prefactor6, invSigSq, shift
  contains
    procedure :: setup => pair_softcore_sf_coul_sf_setup
    procedure :: compute => pair_softcore_sf_coul_sf_compute
    procedure :: mix => pair_softcore_sf_coul_sf_mix
end type pair_softcore_sf_coul_sf

contains

!---------------------------------------------------------------------------------------------------

  subroutine pair_softcore_sf_coul_sf_setup( model, params, iparams )
    class(pair_softcore_sf_coul_sf), intent(inout) :: model
    real(rb), intent(in), optional :: params(:)
    integer,  intent(in), optional :: iparams(:)

    ! Model name:
    model%name = "softcore_sf_coul_sf"

    ! Model parameters:
    model%epsilon = params(1)
    model%sigma = params(2)
    model%lambda = params(3)

    ! Check parameter validity:
    if ((model%lambda < zero).or.(model%lambda > one)) then
      stop "ERROR: invalid parameter lambda in pair_softcore_sf_coul_sf setup"
    end if

    ! Pre-computed quantities:
    model%prefactor = 4.0_rb * model%epsilon * model%lambda**exponent_n
    model%prefactor6 = 6.0_rb * model%prefactor
    model%invSigSq = one/model%sigma**2
    model%shift = alpha*(one - model%lambda)**exponent_p

    ! Activate shifted-force status:
    model%shifted_force_vdw = .true.
    model%shifted_force_coul = .true.

  end subroutine pair_softcore_sf_coul_sf_setup

!---------------------------------------------------------------------------------------------------

  subroutine pair_softcore_sf_coul_sf_compute( model, Eij, Wij, invR2, Qi, Qj )
    class(pair_softcore_sf_coul_sf), intent(in)  :: model
    real(rb),                        intent(out) :: Eij, Wij
    real(rb),                        intent(in)  :: invR2, Qi, Qj

    real(rb) :: rsig2, rsig6, sinv, sinvSq, sinvCb, rFc, invR, QiQj, QiQjbyR

    rsig2 = model%invSigSq/invR2
    rsig6 = rsig2*rsig2*rsig2
    sinv = one/(rsig6 + model%shift)
    sinvSq = sinv*sinv
    sinvCb = sinv*sinvSq
    invR = sqrt(invR2)
    QiQj = Qi*Qj
    QiQjbyR = QiQj*invR
    rFc = (model%fshift_vdw + QiQj*model%fshift_coul)/invR
    Eij = model%prefactor*(sinvSq - sinv) + QiQjbyR + model%eshift_vdw + QiQj*model%eshift_coul + rFc
    Wij = model%prefactor6*rsig6*(sinvCb + sinvCb - sinvSq) + QiQjbyR - rFc

  end subroutine pair_softcore_sf_coul_sf_compute

!---------------------------------------------------------------------------------------------------

  function pair_softcore_sf_coul_sf_mix( this, other ) result( mixed )
    class(pair_softcore_sf_coul_sf), intent(in) :: this
    class(cPairModel),               intent(in) :: other
    class(cPairModel), pointer :: mixed

    select type (other)

      type is (pair_softcore_sf_coul_sf)
        allocate(pair_softcore_sf_coul_sf :: mixed)
        call mixed % setup( [sqrt(this%epsilon*other%epsilon), &
                             half*(this%sigma + other%sigma),  &
                             this%lambda*other%lambda          ] )

!      type is (pair_softcore_sf)
!        allocate(pair_softcore_sf_coul_sf :: mixed)
!        call mixed % setup( [sqrt(this%epsilon*other%epsilon), &
!                             half*(this%sigma + other%sigma),  &
!                             this%lambda*other%lambda          ] )

      type is (pair_lj_sf_coul_sf)
        allocate(pair_softcore_sf_coul_sf :: mixed)
        call mixed % setup( [sqrt(this%epsilon*other%epsilon), &
                             half*(this%sigma + other%sigma),  &
                             this%lambda                       ] )

      type is (pair_lj_sf)
        allocate(pair_softcore_sf_coul_sf :: mixed)
        call mixed % setup( [sqrt(this%epsilon*other%epsilon), &
                             half*(this%sigma + other%sigma),  &
                             this%lambda                       ] )

      class default
        mixed => null()

    end select

  end function pair_softcore_sf_coul_sf_mix

!---------------------------------------------------------------------------------------------------

end module pair_softcore_sf_coul_sf_module
