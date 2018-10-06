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
use nonBondedModelClass

implicit none

!> Abstract class for pair interaction models
type, abstract, extends(cNonBondedModel) :: cPairModel
  contains
    procedure(cPairModel_mix),     deferred :: mix
    procedure :: kspace_setup => cPairModel_kspace_setup
end type cPairModel

abstract interface

  subroutine cPairModel_compute( model, Eij, Wij, invR, invR2 )
    import
    class(cPairModel), intent(in)    :: model
    real(rb),          intent(out)   :: Eij, Wij
    real(rb),          intent(in)    :: invR, invR2
  end subroutine cPairModel_compute

  subroutine cPairModel_energy( model, Eij, invR, invR2 )
    import
    class(cPairModel), intent(in)    :: model
    real(rb),          intent(out)   :: Eij
    real(rb),          intent(in)    :: invR, invR2
  end subroutine cPairModel_energy

  subroutine cPairModel_virial( model, Wij, invR, invR2 )
    import
    class(cPairModel), intent(in)    :: model
    real(rb),          intent(out)   :: Wij
    real(rb),          intent(in)    :: invR, invR2
  end subroutine cPairModel_virial

  function cPairModel_mix( this, other ) result( mixed )
    import
    class(cPairModel), intent(in) :: this, other
    class(cPairModel), pointer    :: mixed
  end function cPairModel_mix

end interface

!> Container structure for pair models
type pairContainer
  class(cPairModel), allocatable :: model
  logical  :: coulomb = .false.
  real(rb) :: kCoul = zero
  contains
    procedure :: assign => pairContainer_assign
    generic :: assignment(=) => assign
    procedure :: mix => pairContainer_mix
end type pairContainer

!> Class definition for pair model "none"
type, extends(cPairModel) :: pair_none
  contains
    procedure :: setup => pair_none_setup
    procedure :: compute => pair_none_compute
    procedure :: energy => pair_none_energy
    procedure :: virial => pair_none_virial
    procedure :: mix => pair_none_mix
end type pair_none

contains

!===================================================================================================
!                             P A I R     M O D E L    C L A S S
!===================================================================================================

  subroutine cPairModel_kspace_setup( model, alpha )
    class(cPairModel), intent(inout) :: model
    real(rb),          intent(in)    :: alpha
  end subroutine cPairModel_kspace_setup

!===================================================================================================
!                         P A I R     M O D E L    C O N T A I N E R
!===================================================================================================

  subroutine pairContainer_assign( new, old )
    class(pairContainer), intent(inout) :: new
    type(modelContainer),      intent(in)    :: old

    if (allocated(new%model)) deallocate( new%model )

    if (allocated(old%model)) then
      select type (model => old%model)
        class is (cPairModel)
          allocate( new%model, source = model )
        class default
          call error( "pair model assignment", trim(model%name)//" is not a pair model" )
      end select
    end if

  end subroutine pairContainer_assign

!---------------------------------------------------------------------------------------------------

  function pairContainer_mix( a, b ) result( c )
    class(pairContainer), intent(in) :: a, b
    type(pairContainer)              :: c

    class(cPairModel), pointer :: mixed

    mixed => b % model % mix( a%model )
    if (.not.associated(mixed)) mixed => a % model % mix( b%model )
    if (associated(mixed)) then
      allocate( c%model, source = mixed )
      deallocate( mixed )
    else
      allocate( c%model, source = pair_none(name="none") )
      call warning( "no mixing rule found for models "//trim(a%model%name)// &
                    " and "//trim(b%model%name) )
    end if

    c%Coulomb = a%Coulomb .and. b%Coulomb
    if (c%Coulomb) c%kCoul = sqrt(a%kCoul*b%kCoul)

  end function pairContainer_mix

!===================================================================================================
!                                   P A I R     N O N E
!===================================================================================================

  type(c_ptr) function EmDee_pair_none() bind(C,name="EmDee_pair_none")
    type(pair_none) :: model
    call model % setup()
    EmDee_pair_none = model % deliver()
  end function EmDee_pair_none

!---------------------------------------------------------------------------------------------------

  subroutine pair_none_setup( model, params, iparams )
    class(pair_none), intent(inout) :: model
    real(rb), intent(in), optional :: params(:)
    integer,  intent(in), optional :: iparams(:)
    model%name = "none"
  end subroutine pair_none_setup

!---------------------------------------------------------------------------------------------------

  subroutine pair_none_compute( model, Eij, Wij, invR, invR2 )
    class(pair_none),  intent(in)  :: model
    real(rb),          intent(out) :: Eij, Wij
    real(rb),          intent(in)  :: invR, invR2
    Eij = zero
    Wij = zero
  end subroutine pair_none_compute

!---------------------------------------------------------------------------------------------------

  subroutine pair_none_energy( model, Eij, invR, invR2 )
    class(pair_none), intent(in)  :: model
    real(rb),         intent(out) :: Eij
    real(rb),         intent(in)  :: invR, invR2
    Eij = zero
  end subroutine pair_none_energy

!---------------------------------------------------------------------------------------------------

  subroutine pair_none_virial( model, Wij, invR, invR2 )
    class(pair_none), intent(in)  :: model
    real(rb),         intent(out) :: Wij
    real(rb),         intent(in)  :: invR, invR2
    Wij = zero
  end subroutine pair_none_virial

!---------------------------------------------------------------------------------------------------

  function pair_none_mix( this, other ) result( mixed )
    class(pair_none),  intent(in) :: this
    class(cPairModel), intent(in) :: other
    class(cPairModel), pointer    :: mixed

    ! Mixing rule: pair_none + any pair model => pair_none
    allocate(pair_none :: mixed)
    call mixed % setup()

  end function pair_none_mix

!===================================================================================================

end module pairModelClass
