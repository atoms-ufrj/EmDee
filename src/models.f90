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
  integer(ib) :: id
  real(rb)    :: p1, p2, p3, p4
  real(rb)    :: f14 = 0.0_rb
end type md_model

type model_ptr
  type(md_model), pointer :: model => null()
end type model_ptr

integer(ib), parameter :: mLJ       = 1, &
                          mLJSF     = 2, &
                          mLJCOULSF = 3, &
                          mHARMOMIC = 4, &
                          mMORSE    = 5

contains

!===================================================================================================
!                                      P A I R     M O D E L S
!===================================================================================================

  type(md_model) function pair_lj( sigma, epsilon ) bind(C)
    real(rb), value :: sigma, epsilon
    pair_lj = md_model( mLJ, sigma*sigma, 4.0_rb*epsilon, 0.0_rb, 0.0_rb )
  end function pair_lj

!---------------------------------------------------------------------------------------------------

  type(md_model) function pair_lj_sf( sigma, epsilon, cutoff ) bind(C)
    real(rb), value :: sigma, epsilon, cutoff
    real(rb) :: sr6, sr12, eps4, Ec, Fc
    sr6 = (sigma/cutoff)**6
    sr12 = sr6*sr6
    eps4 = 4.0_rb*epsilon
    Ec = eps4*(sr12 - sr6)
    Fc = 6.0_rb*(eps4*sr12 + Ec)/cutoff
    pair_lj_sf = md_model( mLJSF, sigma**2, eps4, Fc, -(Ec + Fc*cutoff) )
  end function pair_lj_sf

!---------------------------------------------------------------------------------------------------

  type(md_model) function pair_lj_coul_sf( sigma, epsilon, permittivity, cutoff ) bind(C)
    real(rb), value :: sigma, epsilon, permittivity, cutoff
    real(rb) :: sr6, sr12, eps4, Ec, Fc
    sr6 = (sigma/cutoff)**6
    sr12 = sr6*sr6
    eps4 = 4.0_rb*epsilon
    Ec = eps4*(sr12 - sr6)
    Fc = 6.0_rb*(eps4*sr12 + Ec)/cutoff
    pair_lj_coul_sf = md_model( mLJSF, sigma**2, eps4, Fc, -(Ec + Fc*cutoff) )
  end function pair_lj_coul_sf

!===================================================================================================
!                                      B O N D     M O D E L S
!===================================================================================================

  type(md_model) function bond_harmonic( k, r0 ) bind(C)
    real(rb), value :: k, r0
    bond_harmonic = md_model( mHARMOMIC, r0, -k, 0.5_rb*k, 0.0_rb )
  end function bond_harmonic

!---------------------------------------------------------------------------------------------------

  type(md_model) function bond_morse( D, alpha, r0 ) bind(C)
    real(rb), value :: D, alpha, r0
    bond_morse = md_model( mMORSE, r0, -alpha, D, -2.0_rb*D*alpha )
  end function bond_morse

!===================================================================================================
!                                    A N G L E     M O D E L S
!===================================================================================================

  type(md_model) function angle_harmonic( k, theta0 ) bind(C)
    real(rb), value :: k, theta0
    angle_harmonic = md_model( mHARMOMIC, theta0, -k, 0.5_rb*k, 0.0_rb )
  end function angle_harmonic

!===================================================================================================
!                                 D I H E D R A L     M O D E L S
!===================================================================================================

  type(md_model) function dihedral_harmonic( k, phi0 ) bind(C)
    real(rb), value :: k, phi0
    dihedral_harmonic = md_model( mHARMOMIC, phi0, -k, 0.5_rb*k, 0.0_rb )
  end function dihedral_harmonic

!---------------------------------------------------------------------------------------------------

end module models
