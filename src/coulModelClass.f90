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

module coulModelClass

use global
use modelClass

implicit none

!> Abstract class for coul interaction models
type, abstract, extends(cModel) :: cCoulModel
  logical  :: shifted_force = .false.
  real(rb) :: eshift = zero
  real(rb) :: fshift = zero
  contains
    procedure(cCoulModel_compute), deferred :: compute
    procedure :: shifting_setup => cCoulModel_shifting_setup
end type cCoulModel

abstract interface

  subroutine cCoulModel_compute( model, Eij, Wij, invR2, Qi, Qj )
    import
    class(cCoulModel), intent(in)  :: model
    real(rb),          intent(out) :: Eij, Wij
    real(rb),          intent(in)  :: invR2, Qi, Qj
  end subroutine cCoulModel_compute

end interface

!> Container structure for coul models
type coulModelContainer
  class(cCoulModel), allocatable :: model
  contains
    procedure :: coulModelContainer_assign
    generic :: assignment(=) => coulModelContainer_assign
end type coulModelContainer

!> Class definition for coul model "none"
type, extends(cCoulModel) :: coul_none
  contains
    procedure :: setup => coul_none_setup
    procedure :: compute => coul_none_compute
end type coul_none

contains

!===================================================================================================
!                             C O U L     M O D E L    C L A S S
!===================================================================================================

  subroutine cCoulModel_shifting_setup( model, cutoff )
    class(cCoulModel), intent(inout) :: model
    real(rb),          intent(in)    :: cutoff

    real(rb) :: invR2, E, W

    ! Zero energy and force shifts:
    model%fshift = zero
    model%eshift = zero

    ! Compute energies and virials at cutoff:
    invR2 = one/cutoff**2
    call model%compute( E, W, invR2, one, one )

    ! Update van der Waals energy and force shifts:
    model%fshift = W/cutoff
    model%eshift = -E
    if (model%shifted_force) model%eshift = model%eshift - W

  end subroutine cCoulModel_shifting_setup

!===================================================================================================
!                         C O U L     M O D E L    C O N T A I N E R
!===================================================================================================

  subroutine coulModelContainer_assign( new, old )
    class(coulModelContainer), intent(inout) :: new
    type(modelContainer), intent(in)    :: old

    if (allocated(new%model)) deallocate( new%model )

    if (allocated(old%model)) then
      select type (model => old%model)
        class is (cCoulModel)
          allocate( new%model, source = model )
        class default
          stop "ERROR: cannot assign a coul model type from another model type"
      end select
    end if

  end subroutine coulModelContainer_assign

!===================================================================================================
!                                   C O U L     N O N E
!===================================================================================================

  type(c_ptr) function EmDee_coul_none() bind(C,name="EmDee_coul_none")
    type(coul_none), pointer :: model
    allocate(model)
    call model % setup()
    EmDee_coul_none = model % deliver()
  end function EmDee_coul_none

!---------------------------------------------------------------------------------------------------

  subroutine coul_none_setup( model, params, iparams )
    class(coul_none), intent(inout) :: model
    real(rb), intent(in), optional :: params(:)
    integer,  intent(in), optional :: iparams(:)
    model%name = "none"
  end subroutine coul_none_setup

!---------------------------------------------------------------------------------------------------

  subroutine coul_none_compute( model, Eij, Wij, invR2, Qi, Qj )
    class(coul_none), intent(in)  :: model
    real(rb),         intent(out) :: Eij, Wij
    real(rb),         intent(in)  :: invR2, Qi, Qj
    Eij = zero
    Wij = zero
  end subroutine coul_none_compute

!---------------------------------------------------------------------------------------------------

end module coulModelClass
