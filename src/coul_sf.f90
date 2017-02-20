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

module coul_sf_module

use coulModelClass

implicit none

!> Abstract class for coul model sf
!!
!! NOTES: 1) model parameters must be declared individually and tagged with a comment mark "!<>"
!!        2) recognizable parameter types are real(rb) and integer(ib)
!!        3) allocatable one-dimensional arrays (i.e. vectors) are permitted as parameters
!!        4) an integer(ib) scalar parameter - a size - must necessarily succeed every allocatable
!!           parameter or series of equally-sized allocatable parameters.

type, extends(cCoulModel) :: coul_sf
  contains
    procedure :: setup => coul_sf_setup
    procedure :: compute => coul_sf_compute
    procedure :: energy  => coul_sf_energy
    procedure :: virial  => coul_sf_virial
    procedure :: unsplit => coul_sf_unsplit
end type coul_sf

contains

!---------------------------------------------------------------------------------------------------

  subroutine coul_sf_setup( model, params, iparams )
    class(coul_sf),  intent(inout) :: model
    real(rb), intent(in), optional :: params(:)
    integer,  intent(in), optional :: iparams(:)

    ! Model name:
    model%name = "sf"

    ! Activate shifted-force status:
    model%shifted_force = .true.

  end subroutine coul_sf_setup

!---------------------------------------------------------------------------------------------------

  subroutine coul_sf_compute( model, ECij, WCij, noInvR, invR, invR2, QiQj )
    class(coul_sf), intent(in)    :: model
    real(rb),       intent(out)   :: ECij, WCij
    logical,        intent(inout) :: noInvR
    real(rb),       intent(inout) :: invR
    real(rb),       intent(in)    :: invR2, QiQj

    real(rb) :: rFc, QiQjbyR

    if (noInvR) then
      invR = sqrt(invR2)
      noInvR = .false.
    end if
    QiQjbyR = QiQj*invR
    rFc = QiQj*model%fshift/invR
    ECij = QiQjbyR + QiQj*model%eshift + rFc
    WCij = QiQjbyR - rFc

  end subroutine coul_sf_compute

!---------------------------------------------------------------------------------------------------

  subroutine coul_sf_energy( model, ECij, noInvR, invR, invR2, QiQj )
    class(coul_sf), intent(in)    :: model
    real(rb),       intent(out)   :: ECij
    logical,        intent(inout) :: noInvR
    real(rb),       intent(inout) :: invR
    real(rb),       intent(in)    :: invR2, QiQj

    if (noInvR) then
      invR = sqrt(invR2)
      noInvR = .false.
    end if
    ECij = QiQj*(invR + model%eshift + model%fshift/invR)

  end subroutine coul_sf_energy

!---------------------------------------------------------------------------------------------------

  subroutine coul_sf_virial( model, WCij, noInvR, invR, invR2, QiQj )
    class(coul_sf), intent(in)    :: model
    real(rb),       intent(out)   :: WCij
    real(rb),       intent(inout) :: invR
    logical,        intent(inout) :: noInvR
    real(rb),       intent(in)    :: invR2, QiQj

    if (noInvR) then
      invR = sqrt(invR2)
      noInvR = .false.
    end if
    WCij = QiQj*(invR - model%fshift/invR)

  end subroutine coul_sf_virial

!---------------------------------------------------------------------------------------------------

  subroutine coul_sf_unsplit( model, WCij, noInvR, invR, invR2, QiQj )
    class(coul_sf), intent(in)    :: model
    real(rb),       intent(out)   :: WCij
    real(rb),       intent(inout) :: invR
    logical,        intent(inout) :: noInvR
    real(rb),       intent(in)    :: invR2, QiQj

    if (noInvR) then
      invR = sqrt(invR2)
      noInvR = .false.
    end if
    WCij = QiQj*(invR - model%fshift/invR)

  end subroutine coul_sf_unsplit

!---------------------------------------------------------------------------------------------------

end module coul_sf_module
