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
    procedure :: energy => coul_long_energy
    procedure :: virial => coul_long_virial
end type coul_long

contains

!---------------------------------------------------------------------------------------------------

  subroutine coul_long_setup( model, params, iparams )
    class(coul_long),   intent(inout) :: model
    real(rb), optional, intent(in)    :: params(:)
    integer,  optional, intent(in)    :: iparams(:)

    ! Model name:
    model%name = "long"

    ! Attest that model requires kspace:
    model%requires_kspace = .true.

  end subroutine coul_long_setup

!---------------------------------------------------------------------------------------------------

  subroutine coul_long_kspace_setup( model, alpha )
    class(coul_long), intent(inout) :: model
    real(rb),         intent(in)    :: alpha

    model%alpha = alpha
    model%beta = two*alpha/sqrt(Pi)

  end subroutine coul_long_kspace_setup

!---------------------------------------------------------------------------------------------------

  subroutine coul_long_compute( model, ECij, WCij, noInvR, invR, invR2, QiQj )
    class(coul_long), intent(in)    :: model
    real(rb),         intent(out)   :: ECij, WCij
    logical,          intent(inout) :: noInvR
    real(rb),         intent(inout) :: invR
    real(rb),         intent(in)    :: invR2, QiQj

    real(rb) :: x, expmx2

    if (noInvR) then
      invR = sqrt(invR2)
      noInvR = .false.
    end if
    x = model%alpha/invR
    expmx2 = exp(-x*x)
    ECij = QiQj*uerfc(x,expmx2)*invR
    WCij = ECij + QiQj*model%beta*expmx2

  end subroutine coul_long_compute

!---------------------------------------------------------------------------------------------------

  subroutine coul_long_energy( model, ECij, noInvR, invR, invR2, QiQj )
    class(coul_long), intent(in)    :: model
    real(rb),         intent(out)   :: ECij
    logical,          intent(inout) :: noInvR
    real(rb),         intent(inout) :: invR
    real(rb),         intent(in)    :: invR2, QiQj

    real(rb) :: x

    if (noInvR) then
      invR = sqrt(invR2)
      noInvR = .false.
    end if
    x = model%alpha/invR
    ECij = QiQj*uerfc(x,exp(-x*x))*invR

  end subroutine coul_long_energy

!---------------------------------------------------------------------------------------------------

  subroutine coul_long_virial( model, WCij, noInvR, invR, invR2, QiQj )
    class(coul_long), intent(in)    :: model
    real(rb),         intent(out)   :: WCij
    real(rb),         intent(inout) :: invR
    logical,          intent(inout) :: noInvR
    real(rb),         intent(in)    :: invR2, QiQj

    real(rb) :: x, expmx2

    if (noInvR) then
      invR = sqrt(invR2)
      noInvR = .false.
    end if
    x = model%alpha/invR
    expmx2 = exp(-x*x)
    WCij = QiQj*(uerfc(x,expmx2)*invR + model%beta*expmx2)

  end subroutine coul_long_virial

!---------------------------------------------------------------------------------------------------

end module coul_long_module
