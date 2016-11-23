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

module pairModelClass

use global
use, intrinsic :: iso_c_binding
use modelClass

type, abstract, extends(cModel) :: cPairModel
  contains
    procedure(cPairModel_setup), deferred :: setup
    procedure(cPairModel_compute), deferred :: compute
    procedure(cPairModel_mix), deferred :: mix
end type cPairModel

type, extends(cPairModel) :: pair_none
  contains
    procedure :: setup => pair_none_setup
    procedure :: compute => pair_none_compute
    procedure :: mix => pair_none_mix
end type pair_none

type cPairModelPtr
  class(cPairModel), allocatable :: model
  logical :: overridable
  contains
    procedure :: cPairModelPtr_assign
    generic :: assignment(=) => cPairModelPtr_assign

    procedure :: mix => cPairModelPtr_mix

end type cPairModelPtr

abstract interface

  subroutine cPairModel_setup( model, params )
    import
    class(cPairModel), intent(inout) :: model
    real(rb),          intent(in)    :: params(:)
  end subroutine cPairModel_setup

  subroutine cPairModel_compute( model, Eij, Wij, invR2, Qi, Qj )
    import
    class(cPairModel), intent(in)  :: model
    real(rb),          intent(out) :: Eij, Wij
    real(rb),          intent(in)  :: invR2, Qi, Qj
  end subroutine cPairModel_compute

  function cPairModel_mix( this, other ) result( mixed )
    import
    class(cPairModel), intent(in) :: this, other
    class(cPairModel), pointer    :: mixed
  end function cPairModel_mix

end interface

contains

!---------------------------------------------------------------------------------------------------

!  subroutine cPairModelPtr_assign( new, old )
!    class(cPairModelPtr), intent(inout) :: new
!    type(cPairModelPtr),  intent(in)    :: old

!    if (associated(new%model)) deallocate( new%model )
!    new%model => old%model

!  end subroutine cPairModelPtr_assign

!---------------------------------------------------------------------------------------------------

  type(c_ptr) function EmDee_pair_none() bind(C,name="EmDee_pair_none")
    type(pair_none), pointer :: model
    allocate(model)
    EmDee_pair_none = model % deliver()
  end function EmDee_pair_none

!---------------------------------------------------------------------------------------------------

  subroutine pair_none_setup( model, params )
    class(pair_none), intent(inout) :: model
    real(rb),         intent(in)    :: params(:)
  end subroutine pair_none_setup

!---------------------------------------------------------------------------------------------------

  subroutine pair_none_compute( model, Eij, Wij, invR2, Qi, Qj )
    class(pair_none), intent(in)  :: model
    real(rb),         intent(out) :: Eij, Wij
    real(rb),         intent(in)  :: invR2, Qi, Qj
    Eij = zero
    Wij = zero
  end subroutine pair_none_compute

!---------------------------------------------------------------------------------------------------

  function pair_none_mix( this, other ) result( mixed )
    class(pair_none),  intent(in) :: this
    class(cPairModel), intent(in) :: other
    class(cPairModel), pointer    :: mixed
    allocate(pair_none :: mixed)
  end function pair_none_mix

!---------------------------------------------------------------------------------------------------

  subroutine cPairModelPtr_assign( new, old )
    class(cPairModelPtr), intent(inout) :: new
    type(cModelPtr),      intent(in)    :: old

    if (allocated(new%model)) deallocate( new%model )

    if (allocated(old%model)) then
      select type (model => old%model)
        class is (cPairModel)
          allocate( new%model, source = model )
        class default
          stop "ERROR: cannot assign a pair model type from another model type"
      end select
    end if

  end subroutine cPairModelPtr_assign

!---------------------------------------------------------------------------------------------------

  function cPairModelPtr_mix( a, b ) result( c )
    class(cPairModelPtr), intent(in) :: a
    class(cPairModel),    intent(in) :: b
    type(cPairModelPtr)              :: c
    allocate( pair_none::c%model )
    c%overridable = .true.
  end function cPairModelPtr_mix

!---------------------------------------------------------------------------------------------------

end module pairModelClass
