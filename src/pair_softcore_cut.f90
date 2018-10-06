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

module pair_softcore_cut_module

use global
use pairModelClass
use pair_lj_cut_module

implicit none

real(rb), parameter, private :: alpha      = half, &
                                exponent_n = one,  &
                                exponent_p = one

!> Abstract class for pair model softcore_cut
!!
!! NOTES: 1) model parameters must be declared individually and tagged with a comment mark "!<>"
!!        2) recognizable parameter types are real(rb) and integer(ib)
!!        3) allocatable one-dimensional arrays (i.e. vectors) are permitted as parameters
!!        4) an integer(ib) scalar parameter - a size - must necessarily succeed every allocatable
!!           parameter or series of equally-sized allocatable parameters.

type, extends(cPairModel) :: pair_softcore_cut
  real(rb) :: epsilon !<> Depth of the potential well
  real(rb) :: sigma   !<> Distance at which the potential is zero
  real(rb) :: lambda  !<> Coupling parameter value

  real(rb) :: prefactor, prefactor6, invSigSq, shift
  contains
    procedure :: setup => pair_softcore_cut_setup
    procedure :: compute => pair_softcore_cut_compute
    procedure :: energy => pair_softcore_cut_energy
    procedure :: virial => pair_softcore_cut_virial
    procedure :: mix => pair_softcore_cut_mix
end type pair_softcore_cut

contains

!---------------------------------------------------------------------------------------------------

  subroutine pair_softcore_cut_setup( model, params, iparams )
    class(pair_softcore_cut), intent(inout) :: model
    real(rb), intent(in), optional :: params(:)
    integer,  intent(in), optional :: iparams(:)

    ! Model name:
    model%name = "softcore_cut"

    ! Model parameters:
    model%epsilon = params(1)
    model%sigma = params(2)
    model%lambda = params(3)

    ! Check parameter validity:
    if ((model%lambda < zero).or.(model%lambda > one)) then
      call error( "pair_softcore_cut setup", "out-of-range parameter lambda" )
    end if

    ! Pre-computed quantities:
    model%prefactor = 4.0_rb * model%epsilon * model%lambda**exponent_n
    model%prefactor6 = 6.0_rb * model%prefactor
    model%invSigSq = one/model%sigma**2
    model%shift = alpha*(one - model%lambda)**exponent_p

  end subroutine pair_softcore_cut_setup

!---------------------------------------------------------------------------------------------------

  subroutine pair_softcore_cut_compute( model, Eij, Wij, invR, invR2, QiQj )
    class(pair_softcore_cut), intent(in)  :: model
    real(rb),                 intent(out) :: Eij, Wij
    real(rb),                 intent(in)  :: invR, invR2, QiQj

    real(rb) :: rsig2, rsig6, sinv, sinvSq, sinvCb

    rsig2 = model%invSigSq/invR2
    rsig6 = rsig2*rsig2*rsig2
    sinv = one/(rsig6 + model%shift)
    sinvSq = sinv*sinv
    sinvCb = sinv*sinvSq
    Eij = model%prefactor*(sinvSq - sinv)
    Wij = model%prefactor6*rsig6*(sinvCb + sinvCb - sinvSq)

  end subroutine pair_softcore_cut_compute

!---------------------------------------------------------------------------------------------------

  subroutine pair_softcore_cut_energy( model, Eij, invR, invR2, QiQj )
    class(pair_softcore_cut), intent(in)  :: model
    real(rb),                 intent(out) :: Eij
    real(rb),                 intent(in)  :: invR, invR2, QiQj

    real(rb) :: rsig2, rsig6, sinv

    rsig2 = model%invSigSq/invR2
    rsig6 = rsig2*rsig2*rsig2
    sinv = one/(rsig6 + model%shift)
    Eij = model%prefactor*(sinv*sinv - sinv)

  end subroutine pair_softcore_cut_energy

!---------------------------------------------------------------------------------------------------

  subroutine pair_softcore_cut_virial( model, Wij, invR, invR2, QiQj )
    class(pair_softcore_cut), intent(in)  :: model
    real(rb),                 intent(out) :: Wij
    real(rb),                 intent(in)  :: invR, invR2, QiQj

    real(rb) :: rsig2, rsig6, sinv, sinvSq, sinvCb

    rsig2 = model%invSigSq/invR2
    rsig6 = rsig2*rsig2*rsig2
    sinv = one/(rsig6 + model%shift)
    sinvSq = sinv*sinv
    sinvCb = sinv*sinvSq
    Wij = model%prefactor6*rsig6*(sinvCb + sinvCb - sinvSq)

  end subroutine pair_softcore_cut_virial

!---------------------------------------------------------------------------------------------------

  function pair_softcore_cut_mix( this, other ) result( mixed )
    class(pair_softcore_cut), intent(in) :: this
    class(cPairModel),        intent(in) :: other
    class(cPairModel),        pointer    :: mixed

    select type (other)

      class is (pair_softcore_cut)
        allocate(pair_softcore_cut :: mixed)
        call mixed % setup( [sqrt(this%epsilon*other%epsilon), &
                             half*(this%sigma + other%sigma),  &
                             this%lambda*other%lambda          ] )

      class is (pair_lj_cut)
        allocate(pair_softcore_cut :: mixed)
        call mixed % setup( [sqrt(this%epsilon*other%epsilon), &
                             half*(this%sigma + other%sigma),  &
                             this%lambda                       ] )

      class default
        mixed => null()

    end select

  end function pair_softcore_cut_mix

!---------------------------------------------------------------------------------------------------

end module pair_softcore_cut_module
