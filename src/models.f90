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

type, bind(C) :: tModel
  integer(ib) :: id
  real(rb)    :: p1, p2, p3, p4
  real(rb)    :: f14 = 0.0_rb
end type tModel

type tModelPtr
  type(tModel), pointer :: model => null()
end type tModelPtr

integer(ib), parameter :: mLJ       = 1, &
                          mLJSF     = 2, &
                          mLJCOULSF = 3, &
                          mHARMOMIC = 4, &
                          mMORSE    = 5

contains

!===================================================================================================
!                                      P A I R     M O D E L S
!===================================================================================================

  type(tModel) function pair_lj( sigma, epsilon ) bind(C)
    real(rb), value :: sigma, epsilon
    pair_lj = tModel( mLJ, sigma*sigma, 4.0_rb*epsilon, 0.0_rb, 0.0_rb )
  end function pair_lj

!---------------------------------------------------------------------------------------------------

  type(tModel) function pair_lj_sf( sigma, epsilon, cutoff ) bind(C)
    real(rb), value :: sigma, epsilon, cutoff
    real(rb) :: sr6, sr12, eps4, Ec, Fc
    sr6 = (sigma/cutoff)**6
    sr12 = sr6*sr6
    eps4 = 4.0_rb*epsilon
    Ec = eps4*(sr12 - sr6)
    Fc = 6.0_rb*(eps4*sr12 + Ec)/cutoff
    pair_lj_sf = tModel( mLJSF, sigma**2, eps4, Fc, -(Ec + Fc*cutoff) )
  end function pair_lj_sf

!---------------------------------------------------------------------------------------------------

  type(tModel) function pair_lj_coul_sf( sigma, epsilon, permittivity, cutoff ) bind(C)
    real(rb), value :: sigma, epsilon, permittivity, cutoff
    real(rb) :: sr6, sr12, eps4, Ec, Fc
    sr6 = (sigma/cutoff)**6
    sr12 = sr6*sr6
    eps4 = 4.0_rb*epsilon
    Ec = eps4*(sr12 - sr6)
    Fc = 6.0_rb*(eps4*sr12 + Ec)/cutoff
    pair_lj_coul_sf = tModel( mLJSF, sigma**2, eps4, Fc, -(Ec + Fc*cutoff) )
  end function pair_lj_coul_sf

!===================================================================================================
!                                      B O N D     M O D E L S
!===================================================================================================

  type(tModel) function bond_harmonic( k, r0 ) bind(C)
    real(rb), value :: k, r0
    bond_harmonic = tModel( mHARMOMIC, r0, -k, 0.5_rb*k, 0.0_rb )
  end function bond_harmonic

!---------------------------------------------------------------------------------------------------

  type(tModel) function bond_morse( D, alpha, r0 ) bind(C)
    real(rb), value :: D, alpha, r0
    bond_morse = tModel( mMORSE, r0, -alpha, D, -2.0_rb*D*alpha )
  end function bond_morse

!===================================================================================================
!                                    A N G L E     M O D E L S
!===================================================================================================

  type(tModel) function angle_harmonic( k, theta0 ) bind(C)
    real(rb), value :: k, theta0
    angle_harmonic = tModel( mHARMOMIC, theta0, -k, 0.5_rb*k, 0.0_rb )
  end function angle_harmonic

!===================================================================================================
!                                 D I H E D R A L     M O D E L S
!===================================================================================================

  type(tModel) function dihedral_harmonic( k, phi0 ) bind(C)
    real(rb), value :: k, phi0
    dihedral_harmonic = tModel( mHARMOMIC, phi0, -k, 0.5_rb*k, 0.0_rb )
  end function dihedral_harmonic

!---------------------------------------------------------------------------------------------------

end module models
