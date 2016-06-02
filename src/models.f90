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

module models

use c_binding

implicit none

type, bind(C) :: md_model
  type(c_ptr) :: data
  type(c_ptr) :: params
end type md_model

integer, parameter, private :: sl = 40

type md_data
  character(sl) :: name
  character(sl), allocatable :: param(:)
  real(rb),      allocatable :: value(:)
end type md_data

integer(ib), parameter :: mNONE     = 0, &  ! No model
                          mLJ       = 1, &  ! Lennard-Jones model
                          mLJSF     = 2, &  ! Lennard-Jones (Shifted-Force) model
                          mLJCOULSF = 3, &  ! Lennard-Jones/Coulomb (Shifted-Force) model
                          mHARMOMIC = 4, &  ! Harmonic model
                          mMORSE    = 5     ! Morse model

type md_params
  integer(ib) :: id = mNONE
  real(rb) :: p1 = 0.0_rb
  real(rb) :: p2 = 0.0_rb
  real(rb) :: p3 = 0.0_rb
  real(rb) :: p4 = 0.0_rb
end type md_params

type param_ptr
  type(md_params), pointer :: params => null()
end type param_ptr

type model_ptr
  type(md_model), pointer :: model => null()
end type model_ptr

private :: set_data, set_params

contains

!===================================================================================================
!                                      P A I R     M O D E L S
!===================================================================================================

  type(md_model) function pair_lj( sigma, epsilon ) result( model ) bind(C)
    real(rb), value :: sigma, epsilon

    model%data = set_data( "Lennard-Jones", ["sigma  ","epsilon"], [sigma, epsilon] )
    model%params = set_params( mLJ, sigma*sigma, 4.0_rb*epsilon )

  end function pair_lj

!---------------------------------------------------------------------------------------------------

  type(md_model) function pair_lj_sf( sigma, epsilon, cutoff ) result( model ) bind(C)
    real(rb), value :: sigma, epsilon, cutoff

    real(rb) :: sr6, sr12, eps4, Ec, Fc

    model%data = set_data( "Lennard-Jones (Shifted-Force)", &
                           ["sigma  ","epsilon", "cutoff "], [sigma, epsilon, cutoff] )
    sr6 = (sigma/cutoff)**6
    sr12 = sr6*sr6
    eps4 = 4.0_rb*epsilon
    Ec = eps4*(sr12 - sr6)
    Fc = 6.0_rb*(eps4*sr12 + Ec)/cutoff
    model%params = set_params( mLJSF, sigma**2, eps4, Fc, -(Ec + Fc*cutoff) )

  end function pair_lj_sf

!---------------------------------------------------------------------------------------------------

  type(md_model) function pair_lj_coul_sf( sigma, epsilon, cutoff ) result( model ) bind(C)
    real(rb), value :: sigma, epsilon, cutoff

    real(rb) :: sr6, sr12, eps4, Ec, Fc

    model%data = set_data( "Lennard-Jones/Coulomb (Shifted-Force)", &
                           ["sigma  ","epsilon", "cutoff "], [sigma, epsilon, cutoff] )
    sr6 = (sigma/cutoff)**6
    sr12 = sr6*sr6
    eps4 = 4.0_rb*epsilon
    Ec = eps4*(sr12 - sr6)
    Fc = 6.0_rb*(eps4*sr12 + Ec)/cutoff
    model%params = set_params( mLJSF, sigma**2, eps4, Fc, -(Ec + Fc*cutoff) )

  end function pair_lj_coul_sf

!===================================================================================================
!                                     M I X I N G     R U L E S
!===================================================================================================

  function cross_interaction( imodel, jmodel ) result( ij )
    type(md_model), pointer, intent(in) :: imodel, jmodel
    type(md_model), pointer             :: ij

    type(md_data),   pointer :: idata, jdata
    type(md_params), pointer :: iparams, jparams

    ij => null()
    if (associated(imodel).and.associated(jmodel)) then
      call c_f_pointer( imodel%data, idata )
      call c_f_pointer( jmodel%data, jdata )
      call c_f_pointer( imodel%params, iparams )
      call c_f_pointer( jmodel%params, jparams )
      if (iparams%id == jparams%id) then
        select case (iparams%id)
          case (mLJ)
            allocate( ij )
            ij = pair_lj( arithmetic(1), geometric(2) )
          case (mLJSF)
            allocate( ij )
            ij = pair_lj_sf( arithmetic(1), geometric(2), arithmetic(3) )
          case (mLJCOULSF)
            allocate( ij )
            ij = pair_lj_coul_sf( arithmetic(1), geometric(2), arithmetic(3) )
        end select
      end if
    end if

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      real(rb) function arithmetic( k )
        integer(ib), intent(in) :: k
        arithmetic = 0.5_rb*(idata%value(k) + jdata%value(k))
      end function arithmetic
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      real(rb) function geometric( k )
        integer(ib), intent(in) :: k
        geometric = sqrt(idata%value(k)*jdata%value(k))
      end function geometric
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end function cross_interaction

!===================================================================================================
!                                      B O N D     M O D E L S
!===================================================================================================

  type(md_model) function bond_harmonic( k, r0 ) result( model ) bind(C)
    real(rb), value :: k, r0

    model%data = set_data( "Harmonic", ["k ", "r0"], [k, r0] )
    model%params = set_params( mHARMOMIC, r0, -k, 0.5_rb*k )

  end function bond_harmonic

!---------------------------------------------------------------------------------------------------

  type(md_model) function bond_morse( D, alpha, r0 ) result( model ) bind(C)
    real(rb), value :: D, alpha, r0

    model%data = set_data( "Morse", ["D    ", "alpha", "r0   "], [D, alpha, r0] )
    model%params = set_params( mMORSE, r0, -alpha, D, -2.0_rb*D*alpha )

  end function bond_morse

!===================================================================================================
!                                    A N G L E     M O D E L S
!===================================================================================================

  type(md_model) function angle_harmonic( k, theta0 ) result( model ) bind(C)
    real(rb), value :: k, theta0

    model%data = set_data( "Harmonic", ["k     ", "theta0"], [k, theta0] )
    model%params = set_params( mHARMOMIC, theta0, -k, 0.5_rb*k )

  end function angle_harmonic

!===================================================================================================
!                                 D I H E D R A L     M O D E L S
!===================================================================================================

  type(md_model) function dihedral_harmonic( k, phi0 ) result( model ) bind(C)
    real(rb), value :: k, phi0

    model%data = set_data( "Harmonic", ["k   ", "phi0"], [k, phi0] )
    model%params = set_params( mHARMOMIC, 0.0_rb, phi0, -k, 0.5_rb*k )

  end function dihedral_harmonic

!===================================================================================================
!                             A U X I L I A R Y     P R O C E D U R E S
!===================================================================================================

  function set_data( name, params, values ) result( data )
    character(*), intent(in) :: name, params(:)
    real(rb),     intent(in) :: values(:)
    type(c_ptr)              :: data

    type(md_data), pointer :: ptr

    allocate( ptr )
    ptr%name = name
    ptr%param = params
    ptr%value = values
    data = c_loc(ptr)

  end function set_data

!---------------------------------------------------------------------------------------------------

  function set_params( id, p1, p2, p3, p4 ) result( params )
    integer(ib), intent(in)           :: id
    real(rb),    intent(in), optional :: p1, p2, p3, p4
    type(c_ptr)                       :: params

    type(md_params), pointer :: ptr

    allocate( ptr )
    ptr%id = id
    if (present(p1)) ptr%p1 = p1
    if (present(p2)) ptr%p2 = p2
    if (present(p3)) ptr%p3 = p3
    if (present(p4)) ptr%p4 = p4
    params = c_loc(ptr)

  end function set_params

!===================================================================================================

end module models
