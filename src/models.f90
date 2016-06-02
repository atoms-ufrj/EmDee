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
  real(rb)    :: p1 = 0.0_rb
  real(rb)    :: p2 = 0.0_rb
  real(rb)    :: p3 = 0.0_rb
  real(rb)    :: p4 = 0.0_rb
end type md_params

type param_ptr
  type(md_params), pointer :: model => null()
end type param_ptr

type data_ptr
  type(md_data), pointer :: data => null()
end type data_ptr

contains

!===================================================================================================
!                                      P A I R     M O D E L S
!===================================================================================================

  type(md_model) function pair_lj( sigma, epsilon ) bind(C)
    real(rb), value :: sigma, epsilon
    pair_lj%data = set_data( "Lennard-Jones", ["sigma  ","epsilon"], [sigma, epsilon] )
    pair_lj%params = set_params( mLJ, sigma*sigma, 4.0_rb*epsilon )
  end function pair_lj

!---------------------------------------------------------------------------------------------------

  type(md_model) function pair_lj_sf( sigma, epsilon, cutoff ) bind(C)
    real(rb), value :: sigma, epsilon, cutoff
    real(rb) :: sr6, sr12, eps4, Ec, Fc
    pair_lj_sf%data = set_data( "Lennard-Jones (Shifted-Force)", &
                                ["sigma  ","epsilon", "cutoff "], [sigma, epsilon, cutoff] )
    sr6 = (sigma/cutoff)**6
    sr12 = sr6*sr6
    eps4 = 4.0_rb*epsilon
    Ec = eps4*(sr12 - sr6)
    Fc = 6.0_rb*(eps4*sr12 + Ec)/cutoff
    pair_lj_sf%params = set_params( mLJSF, sigma**2, eps4, Fc, -(Ec + Fc*cutoff) )
  end function pair_lj_sf

!---------------------------------------------------------------------------------------------------

  type(md_model) function pair_lj_coul_sf( sigma, epsilon, cutoff ) bind(C)
    real(rb), value :: sigma, epsilon, cutoff
    real(rb) :: sr6, sr12, eps4, Ec, Fc
    pair_lj_coul_sf%data = set_data( "Lennard-Jones/Coulomb (Shifted-Force)", &
                                     ["sigma  ","epsilon", "cutoff "], [sigma, epsilon, cutoff] )
    sr6 = (sigma/cutoff)**6
    sr12 = sr6*sr6
    eps4 = 4.0_rb*epsilon
    Ec = eps4*(sr12 - sr6)
    Fc = 6.0_rb*(eps4*sr12 + Ec)/cutoff
    pair_lj_coul_sf%params = set_params( mLJSF, sigma**2, eps4, Fc, -(Ec + Fc*cutoff) )
  end function pair_lj_coul_sf

!===================================================================================================
!                                      B O N D     M O D E L S
!===================================================================================================

  type(md_model) function bond_harmonic( k, r0 ) bind(C)
    real(rb), value :: k, r0
    bond_harmonic%data = set_data( "Harmonic", ["k ", "r0"], [k, r0] )
    bond_harmonic%params = set_params( mHARMOMIC, r0, -k, 0.5_rb*k )
  end function bond_harmonic

!---------------------------------------------------------------------------------------------------

  type(md_model) function bond_morse( D, alpha, r0 ) bind(C)
    real(rb), value :: D, alpha, r0
    bond_morse%data = set_data( "Morse", ["D    ", "alpha", "r0   "], [D, alpha, r0] )
    bond_morse%params = set_params( mMORSE, r0, -alpha, D, -2.0_rb*D*alpha )
  end function bond_morse

!===================================================================================================
!                                    A N G L E     M O D E L S
!===================================================================================================

  type(md_model) function angle_harmonic( k, theta0 ) bind(C)
    real(rb), value :: k, theta0
    angle_harmonic%data = set_data( "Harmonic", ["k     ", "theta0"], [k, theta0] )
    angle_harmonic%params = set_params( mHARMOMIC, theta0, -k, 0.5_rb*k )
  end function angle_harmonic

!===================================================================================================
!                                 D I H E D R A L     M O D E L S
!===================================================================================================

  type(md_model) function dihedral_harmonic( k, phi0 ) bind(C)
    real(rb), value :: k, phi0
    dihedral_harmonic%data = set_data( "Harmonic", ["k   ", "phi0"], [k, phi0] )
    dihedral_harmonic%params = set_params( mHARMOMIC, 0.0_rb, phi0, -k, 0.5_rb*k )
  end function dihedral_harmonic

!===================================================================================================
!                             A U X I L I A R Y     P R O C E D U R E S
!===================================================================================================

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

!---------------------------------------------------------------------------------------------------

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

!===================================================================================================

end module models
