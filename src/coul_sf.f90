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

!> Abstract class for coul model sf
!! NOTES: 1) model parameters must be declared individually and tagged with a comment mark "!<>"
!!        2) allocatable parameters are permitted only for rank=1
!!        3) a series of rank-1 allocatable parameters must be succeeded by an integer parameter,
!!           which will contain their (common) size after allocation
type, extends(cCoulModel) :: coul_sf
  contains
    procedure :: setup => coul_sf_setup
    procedure :: compute => coul_sf_compute
end type coul_sf

contains

!---------------------------------------------------------------------------------------------------

  subroutine coul_sf_setup( model, params, iparams )
    class(coul_sf), intent(inout) :: model
    real(rb), intent(in), optional :: params(:)
    integer,  intent(in), optional :: iparams(:)

    ! Model name:
    model%name = "sf"

    ! Activate shifted-force status:
    model%shifted_force = .true.

  end subroutine coul_sf_setup

!---------------------------------------------------------------------------------------------------

  subroutine coul_sf_compute( model, Eij, Wij, invR2, Qi, Qj )
    class(coul_sf), intent(in)  :: model
    real(rb),       intent(out) :: Eij, Wij
    real(rb),       intent(in)  :: invR2, Qi, Qj

    real(rb) :: rFc, invR, QiQj, QiQjbyR

    invR = sqrt(invR2)
    QiQj = Qi*Qj
    QiQjbyR = QiQj*invR
    rFc = QiQj*model%fshift/invR
    Eij = QiQjbyR + QiQj*model%eshift + rFc
    Wij = QiQjbyR - rFc

  end subroutine coul_sf_compute

!---------------------------------------------------------------------------------------------------

end module coul_sf_module
