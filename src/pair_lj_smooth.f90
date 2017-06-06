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

module pair_lj_smooth_module

use global
use pairModelClass

implicit none

!> Abstract class for pair model lj_smooth
!!
!! NOTES: 1) model parameters must be declared individually and tagged with a comment mark "!<>"
!!        2) recognizable parameter types are real(rb) and integer(ib)
!!        3) allocatable one-dimensional arrays (i.e. vectors) are permitted as parameters
!!        4) an integer(ib) scalar parameter - a size - must necessarily succeed every allocatable
!!           parameter or series of equally-sized allocatable parameters.

type, extends(cPairModel) :: pair_lj_smooth
  real(rb) :: epsilon !<> Depth of the potential well
  real(rb) :: sigma   !<> Distance at which the potential is zero
  real(rb) :: Rm      !<> Distance at which smoothing starts
  real(rb) :: Rc      !<> Cut-off distance

  real(rb) :: eps4, eps24, sigsq, Rm2, invRm, factor
  contains
    procedure :: setup => pair_lj_smooth_setup
    procedure :: compute => pair_lj_smooth_compute
    procedure :: energy => pair_lj_smooth_energy
    procedure :: virial => pair_lj_smooth_virial
    procedure :: mix => pair_lj_smooth_mix
end type pair_lj_smooth

contains

!---------------------------------------------------------------------------------------------------

  subroutine pair_lj_smooth_setup( model, params, iparams )
    class(pair_lj_smooth), intent(inout) :: model
    real(rb), intent(in), optional :: params(:)
    integer,  intent(in), optional :: iparams(:)

    ! Model name:
    model%name = "lj_smooth"

    ! Model parameters:
    model%epsilon = params(1)
    model%sigma = params(2)
    model%Rm = params(3)
    model%Rc = params(4)

    ! Pre-computed quantities:
    model%eps4  = 4.0_rb*model%epsilon
    model%eps24 = 24.0_rb*model%epsilon
    model%sigsq = model%sigma**2
    model%Rm2 = model%Rm**2
    model%invRm = one/model%Rm
    model%factor = one/(model%Rc**2 - model%Rm2)

  end subroutine pair_lj_smooth_setup

!---------------------------------------------------------------------------------------------------

  subroutine pair_lj_smooth_compute( model, Eij, Wij, invR, invR2 )
    class(pair_lj_smooth), intent(in)  :: model
    real(rb),           intent(out) :: Eij, Wij
    real(rb),           intent(in)  :: invR, invR2

    real(rb) :: sr2, sr6, sr12, r2, u, u2, u3, G, WG

    sr2 = model%sigSq*invR2
    sr6 = sr2*sr2*sr2
    sr12 = sr6*sr6
    Eij = model%eps4*(sr12 - sr6)
    Wij = model%eps24*(sr12 + sr12 - sr6)

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

  end subroutine pair_lj_smooth_compute

!---------------------------------------------------------------------------------------------------

  subroutine pair_lj_smooth_energy( model, Eij, invR, invR2 )
    class(pair_lj_smooth), intent(in)  :: model
    real(rb),           intent(out) :: Eij
    real(rb),           intent(in)  :: invR, invR2

    real(rb) :: sr2, sr6, sr12, r2, u, u2, u3, G

    sr2 = model%sigSq*invR2
    sr6 = sr2*sr2*sr2
    sr12 = sr6*sr6
    Eij = model%eps4*(sr12 - sr6)

    if (invR < model%invRm) then
      r2 = one/invR2
      u = model%factor*(r2 - model%Rm2)
      u2 = u*u
      u3 = u*u2
      G = 1.0_rb + u3*(15.0_rb*u - 6.0_rb*u2 - 10.0_rb)
      Eij = Eij*G
    end if

  end subroutine pair_lj_smooth_energy

!---------------------------------------------------------------------------------------------------

  subroutine pair_lj_smooth_virial( model, Wij, invR, invR2 )
    class(pair_lj_smooth), intent(in)  :: model
    real(rb),           intent(out) :: Wij
    real(rb),           intent(in)  :: invR, invR2

    real(rb) :: sr2, sr6, sr12, r2, u, u2, u3, G, WG

    sr2 = model%sigSq*invR2
    sr6 = sr2*sr2*sr2
    sr12 = sr6*sr6
    Wij = model%eps24*(sr12 + sr12 - sr6)

    if (invR < model%invRm) then
      r2 = one/invR2
      u = model%factor*(r2 - model%Rm2)
      u2 = u*u
      u3 = u*u2
      G = 1.0_rb + u3*(15.0_rb*u - 6.0_rb*u2 - 10.0_rb)
      WG = -60.0_rb*u2*(2.0_rb*u - u2 - 1.0_rb)*model%factor*r2
      Wij = Wij*G + model%eps4*(sr12 - sr6)*WG
    end if

  end subroutine pair_lj_smooth_virial

!---------------------------------------------------------------------------------------------------

  function pair_lj_smooth_mix( this, other ) result( mixed )
    class(pair_lj_smooth), intent(in) :: this
    class(cPairModel),  intent(in) :: other
    class(cPairModel),  pointer    :: mixed

    select type (other)

      class is (pair_lj_smooth)
        allocate(pair_lj_smooth :: mixed)
        call mixed % setup( [sqrt(this%epsilon*other%epsilon), &
                             half*(this%sigma + other%sigma),  &
                             half*(this%Rm + other%Rm),        &
                             half*(this%Rc + other%Rc)]        )
      class default
        mixed => null()

    end select

  end function pair_lj_smooth_mix

!---------------------------------------------------------------------------------------------------

end module pair_lj_smooth_module
