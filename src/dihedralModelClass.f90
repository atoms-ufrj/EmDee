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

module dihedralModelClass

use global
use modelClass

implicit none

!> An abstract class for dihedral interaction models:
type, abstract, extends(cModel) :: cDihedralModel
  real(rb) :: factor14 = zero
  contains
    procedure(cDihedralModel_compute), deferred :: compute
end type cDihedralModel

!> A class for no-interaction dihedral model:
type, extends(cDihedralModel) :: dihedral_none
  contains
    procedure :: setup => dihedral_none_setup
    procedure :: compute => dihedral_none_compute
end type dihedral_none

abstract interface

  subroutine cDihedralModel_compute( model, Ed, Fd, phi )
    import
    class(cDihedralModel), intent(in)  :: model
    real(rb),              intent(out) :: Ed, Fd
    real(rb),              intent(in)  :: phi
  end subroutine cDihedralModel_compute

end interface

contains

!===================================================================================================
!                                  D I H E D R A L     N O N E
!===================================================================================================

  type(c_ptr) function EmDee_dihedral_none() bind(C,name="EmDee_dihedral_none")
    type(dihedral_none), pointer :: model
    allocate(model)
    call model % setup()
    EmDee_dihedral_none = model % deliver()
  end function EmDee_dihedral_none

!---------------------------------------------------------------------------------------------------

  subroutine dihedral_none_setup( model, params, iparams )
    class(dihedral_none), intent(inout) :: model
    real(rb), intent(in), optional :: params(:)
    integer,  intent(in), optional :: iparams(:)
    model%name = "none"
  end subroutine dihedral_none_setup

!---------------------------------------------------------------------------------------------------

  subroutine dihedral_none_compute( model, Ed, Fd, phi )
    class(dihedral_none), intent(in)  :: model
    real(rb),             intent(out) :: Ed, Fd
    real(rb),             intent(in)  :: phi
    Ed = zero
    Fd = zero
  end subroutine dihedral_none_compute

!---------------------------------------------------------------------------------------------------

end module dihedralModelClass
