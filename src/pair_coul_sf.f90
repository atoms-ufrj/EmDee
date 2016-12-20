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

module pair_coul_sf_module

use pairModelClass

!> Abstract class for pair model coul_sf
!! NOTES: 1) model parameters must be declared individually and tagged with a comment mark "!<>"
!!        2) allocatable parameters are permitted only for rank=1
!!        3) a series of rank-1 allocatable parameters must be succeeded by an integer parameter,
!!           which will contain their (common) size after allocation
type, extends(cPairModel) :: pair_coul_sf
!  integer :: N !<>
!  integer, allocatable :: A(:) !<>
!  real, allocatable :: B(:) !<>
!  integer :: M !<>
!  real, allocatable :: C(:) !<>
!  real :: D !<>
  contains
    procedure :: setup => pair_coul_sf_setup
    procedure :: compute => pair_coul_sf_compute
    procedure :: mix => pair_coul_sf_mix
end type pair_coul_sf

contains

!---------------------------------------------------------------------------------------------------

  subroutine pair_coul_sf_setup( model, params, iparams )
    class(pair_coul_sf), intent(inout) :: model
    real(rb), intent(in), optional :: params(:)
    integer,  intent(in), optional :: iparams(:)

    ! Model name:
    model%name = "coul_sf"

    ! Activate shifted-force status:
    model%shifted_force_coul = .true.

  end subroutine pair_coul_sf_setup

!---------------------------------------------------------------------------------------------------

  subroutine pair_coul_sf_compute( model, Eij, Wij, invR2, Qi, Qj )
    class(pair_coul_sf), intent(in)  :: model
    real(rb),       intent(out) :: Eij, Wij
    real(rb),       intent(in)  :: invR2, Qi, Qj

    real(rb) :: rFc, invR, QiQj, QiQjbyR

    invR = sqrt(invR2)
    QiQj = Qi*Qj
    QiQjbyR = QiQj*invR
    rFc = QiQj*model%fshift_coul/invR
    Eij = QiQjbyR + QiQj*model%eshift_coul + rFc
    Wij = QiQjbyR - rFc

  end subroutine pair_coul_sf_compute

!---------------------------------------------------------------------------------------------------

  function pair_coul_sf_mix( this, other ) result( mixed )
    class(pair_coul_sf),    intent(in) :: this
    class(cPairModel), intent(in) :: other
    class(cPairModel), pointer :: mixed

    allocate(pair_coul_sf :: mixed)
    call mixed % setup( )

  end function pair_coul_sf_mix

!---------------------------------------------------------------------------------------------------

end module pair_coul_sf_module
