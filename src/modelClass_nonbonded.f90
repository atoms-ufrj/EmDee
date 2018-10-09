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

module nonBondedModelClass

use global
use modelClass

implicit none

integer, parameter :: NONE = 0,                    &
                      SHIFTED = 1,                 &
                      SHIFTED_FORCE = 2,           &
                      SMOOTHED = 3,                &
                      SHIFTED_SMOOTHED = 4,        &
                      SQUARE_SMOOTHED = 5,         &
                      SHIFTED_SQUARE_SMOOTHED = 6

!> Abstract class for non-bonded interaction models
type, abstract, extends(cModel) :: cNonBondedModel
  integer  :: modifier = NONE
  real(rb) :: eshift = zero
  real(rb) :: fshift = zero
  real(rb) :: skin = zero
  real(rb) :: Rm = zero
  real(rb) :: RmSq = zero
  real(rb) :: factor = zero
  real(rb) :: Rm2fac = zero
  contains
    procedure(cNonBondedModel_compute), deferred :: compute
    procedure(cNonBondedModel_energy),  deferred :: energy
    procedure(cNonBondedModel_virial),  deferred :: virial
    procedure :: modifier_setup => cNonBondedModel_modifier_setup
end type cNonBondedModel

abstract interface

  subroutine cNonBondedModel_compute( model, Eij, Wij, invR, invR2 )
    import
    class(cNonBondedModel), intent(in)  :: model
    real(rb),          intent(out) :: Eij, Wij
    real(rb),          intent(in)  :: invR, invR2
  end subroutine cNonBondedModel_compute

  subroutine cNonBondedModel_energy( model, Eij, invR, invR2 )
    import
    class(cNonBondedModel), intent(in)  :: model
    real(rb),          intent(out) :: Eij
    real(rb),          intent(in)  :: invR, invR2
  end subroutine cNonBondedModel_energy

  subroutine cNonBondedModel_virial( model, Wij, invR, invR2 )
    import
    class(cNonBondedModel), intent(in)  :: model
    real(rb),          intent(out) :: Wij
    real(rb),          intent(in)  :: invR, invR2
  end subroutine cNonBondedModel_virial

end interface

contains

!===================================================================================================
!                           N O N B O N D E D     M O D E L    C L A S S
!===================================================================================================

  type(c_ptr) function EmDee_shifted( model ) &
    bind(C,name="EmDee_shifted")
    type(c_ptr), value :: model

    type(modelContainer),  pointer :: container
    class(cNonBondedModel), allocatable :: new

    if (.not.c_associated(model)) then
      call error( "shifted potential assignment", "a valid pair model must be provided" )
    end if
    call c_f_pointer( model, container )

    select type ( pair => container%model )
      class is (cNonBondedModel)
        allocate( new, source = pair )
        new%modifier = SHIFTED
        EmDee_shifted = new % deliver()
      class default
        call error( "shifted potential assignment", "a valid pair model must be provided" )
    end select

  end function EmDee_shifted

!---------------------------------------------------------------------------------------------------

  type(c_ptr) function EmDee_shifted_force( model ) &
    bind(C,name="EmDee_shifted_force")
    type(c_ptr), value :: model

    type(modelContainer),  pointer :: container
    class(cNonBondedModel), allocatable :: new

    if (.not.c_associated(model)) then
      call error( "shifted-force potential assignment", "a valid pair model must be provided" )
    end if
    call c_f_pointer( model, container )

    select type ( pair => container%model )
      class is (cNonBondedModel)
        allocate( new, source = pair )
        new%modifier = SHIFTED_FORCE
        EmDee_shifted_force = new % deliver()
      class default
        call error( "shifted-force potential assignment", "a valid pair model must be provided" )
    end select

  end function EmDee_shifted_force

!---------------------------------------------------------------------------------------------------

  type(c_ptr) function EmDee_smoothed( model, skin ) &
    bind(C,name="EmDee_smoothed")
    type(c_ptr), value :: model
    real(rb),    value :: skin

    type(modelContainer),  pointer :: container
    class(cNonBondedModel), allocatable :: new

    if (.not.c_associated(model)) then
      call error( "smoothed potential assignment", "a valid pair model must be provided" )
    end if
    call c_f_pointer( model, container )

    select type ( pair => container%model )
      class is (cNonBondedModel)
        allocate( new, source = pair )
        new%modifier = SMOOTHED
        new%skin = skin
        EmDee_smoothed = new % deliver()
      class default
        call error( "smoothed potential assignment", "a valid pair model must be provided" )
    end select

  end function EmDee_smoothed

!---------------------------------------------------------------------------------------------------

  type(c_ptr) function EmDee_shifted_smoothed( model, skin ) &
    bind(C,name="EmDee_shifted_smoothed")
    type(c_ptr), value :: model
    real(rb),    value :: skin

    type(modelContainer),  pointer :: container
    class(cNonBondedModel), allocatable :: new

    if (.not.c_associated(model)) then
      call error( "shifted-smoothed potential assignment", &
                  "a valid pair model must be provided" )
    end if
    call c_f_pointer( model, container )

    select type ( pair => container%model )
      class is (cNonBondedModel)
        allocate( new, source = pair )
        new%modifier = SHIFTED_SMOOTHED
        new%skin = skin
        EmDee_shifted_smoothed = new % deliver()
      class default
        call error( "shifted-smoothed potential assignment", &
                    "a valid pair model must be provided" )
    end select

  end function EmDee_shifted_smoothed

!---------------------------------------------------------------------------------------------------

  type(c_ptr) function EmDee_square_smoothed( model, skin ) &
    bind(C,name="EmDee_square_smoothed")
    type(c_ptr), value :: model
    real(rb),    value :: skin

    type(modelContainer),  pointer :: container
    class(cNonBondedModel), allocatable :: new

    if (.not.c_associated(model)) then
      call error( "square-smoothed potential assignment", "a valid pair model must be provided" )
    end if
    call c_f_pointer( model, container )

    select type ( pair => container%model )
      class is (cNonBondedModel)
        allocate( new, source = pair )
        new%modifier = SQUARE_SMOOTHED
        new%skin = skin
        EmDee_square_smoothed = new % deliver()
      class default
        call error( "square-smoothed potential assignment", "a valid pair model must be provided" )
    end select

  end function EmDee_square_smoothed

!---------------------------------------------------------------------------------------------------

  type(c_ptr) function EmDee_shifted_square_smoothed( model, skin ) &
    bind(C,name="EmDee_shifted_square_smoothed")
    type(c_ptr), value :: model
    real(rb),    value :: skin

    type(modelContainer),  pointer :: container
    class(cNonBondedModel), allocatable :: new

    if (.not.c_associated(model)) then
      call error( "shifted-square-smoothed potential assignment", &
                  "a valid pair model must be provided" )
    end if
    call c_f_pointer( model, container )

    select type ( pair => container%model )
      class is (cNonBondedModel)
        allocate( new, source = pair )
        new%modifier = SHIFTED_SQUARE_SMOOTHED
        new%skin = skin
        EmDee_shifted_square_smoothed = new % deliver()
      class default
        call error( "shifted-square-smoothed potential assignment", &
                    "a valid pair model must be provided" )
    end select

  end function EmDee_shifted_square_smoothed

!---------------------------------------------------------------------------------------------------

  subroutine cNonBondedModel_modifier_setup( model, cutoff )
    class(cNonBondedModel), intent(inout) :: model
    real(rb),               intent(in)    :: cutoff

    real(rb) :: invR2, E, W, invR

    ! Zero energy and force shifts:
    model%fshift = zero
    model%eshift = zero

    ! Store cutoff-dependent parameters:
    model%Rm = cutoff - model%skin
    model%RmSq = model%Rm**2

    select case (model%modifier)
      case (SHIFTED, SHIFTED_FORCE, SHIFTED_SMOOTHED, SHIFTED_SQUARE_SMOOTHED)

        ! Compute energies and virials at cutoff:
        invR = one/cutoff
        invR2 = invR*invR
        call model%compute( E, W, invR, invR2 )

        ! Update energy and force shifts:
        if (model%modifier == SHIFTED_FORCE) then
          model%eshift = -(E + W)
          model%fshift = W/cutoff
        else
          model%fshift = zero
          model%eshift = -E
          if (model%modifier == SHIFTED_SMOOTHED) then
            model%factor = one/(cutoff - model%Rm)
            model%Rm2fac = model%factor*model%Rm
          else if (model%modifier == SHIFTED_SQUARE_SMOOTHED) then
            model%factor = one/(cutoff**2 - model%RmSq)
            model%Rm2fac = model%factor*model%RmSq
          end if
        end if

      case (SMOOTHED)

        model%factor = one/(cutoff - model%Rm)
        model%Rm2fac = model%factor*model%Rm

      case (SQUARE_SMOOTHED)

        model%factor = one/(cutoff**2 - model%RmSq)
        model%Rm2fac = model%factor*model%RmSq

    end select

  end subroutine cNonBondedModel_modifier_setup

!---------------------------------------------------------------------------------------------------

end module nonBondedModelClass
