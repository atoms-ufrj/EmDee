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

integer(ib), parameter :: mCOULOMB = 10000

integer(ib), parameter :: mNONE     = 0, &  ! No model
                          mLJ       = 1, &  ! Lennard-Jones
                          mLJSF     = 2, &  ! Lennard-Jones, Shifted-Force
                          mHARMOMIC = 3, &  ! Harmonic
                          mMORSE    = 4     ! Morse

real(rb), parameter, private :: Deg2Rad = 3.14159265358979324_rb/180_rb

type, bind(C) :: EmDee_Model
  integer(ib) :: id = mNONE
  type(c_ptr) :: data = c_null_ptr
  real(rb)    :: p1 = 0.0_rb
  real(rb)    :: p2 = 0.0_rb
  real(rb)    :: p3 = 0.0_rb
  real(rb)    :: p4 = 0.0_rb
  integer(ib) :: external = 1
end type EmDee_Model

type model_ptr
  type(EmDee_Model), pointer :: model => null()
end type model_ptr

private :: set_data

contains

!===================================================================================================
!                                      P A I R     M O D E L S
!===================================================================================================

  function EmDee_pair_none() result(model) bind(C,name="EmDee_pair_none")
    type(EmDee_Model) :: model
    model = EmDee_Model( mNONE )
  end function EmDee_pair_none

!---------------------------------------------------------------------------------------------------

  function EmDee_pair_lj( epsilon, sigma ) result(model) bind(C,name="EmDee_pair_lj")
    real(rb), value   :: epsilon, sigma
    type(EmDee_Model) :: model
    model = EmDee_Model( mLJ, set_data( [epsilon, sigma] ), sigma*sigma, 4.0_rb*epsilon )
  end function EmDee_pair_lj

!---------------------------------------------------------------------------------------------------

  function EmDee_pair_lj_sf( epsilon, sigma, cutoff ) result(model) bind(C,name="EmDee_pair_lj_sf")
    real(rb), value   :: epsilon, sigma, cutoff
    type(EmDee_Model) :: model
    real(rb) :: sr6, sr12, eps4, Ec, Fc
    sr6 = (sigma/cutoff)**6
    sr12 = sr6*sr6
    eps4 = 4.0_rb*epsilon
    Ec = eps4*(sr12 - sr6)
    Fc = 6.0_rb*(eps4*sr12 + Ec)/cutoff
    Ec = -(Ec + Fc*cutoff)
    model = EmDee_Model( mLJSF, set_data( [epsilon, sigma, cutoff] ), sigma**2, eps4, Fc, Ec )
  end function EmDee_pair_lj_sf

!===================================================================================================
!                                     M I X I N G     R U L E S
!===================================================================================================

  function cross_pair( imodel, jmodel ) result( ij )
    type(EmDee_Model), pointer, intent(in) :: imodel, jmodel
    type(EmDee_Model), pointer             :: ij

    real(rb), pointer :: idata(:), jdata(:)

    if (associated(imodel).and.associated(jmodel)) then
      allocate( ij )
      if (match(mLJ,mLJ)) then
        ij = EmDee_pair_lj( geometric(1), arithmetic(2) )
      else if (match(mLJSF,mLJSF)) then
        ij = EmDee_pair_lj_sf( geometric(1), arithmetic(2), arithmetic(3) )
      else if (match(mLJ,mNONE).or.match(mLJSF,mNONE)) then
        ij = EmDee_pair_none()
      else
        deallocate( ij )
        ij => null()
      end if
    else
      ij => null()
    end if

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      logical function match( a, b )
        integer(ib), intent(in) :: a, b
        match = ((imodel%id == a).and.(jmodel%id == b)).or.((imodel%id == b).and.(jmodel%id == a))
      end function match
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      real(rb) function arithmetic( k )
        integer(ib), intent(in) :: k
        call c_f_pointer( imodel%data, idata, [k] )
        call c_f_pointer( jmodel%data, jdata, [k] )
        arithmetic = 0.5_rb*(idata(k) + jdata(k))
      end function arithmetic
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      real(rb) function geometric( k )
        integer(ib), intent(in) :: k
        call c_f_pointer( imodel%data, idata, [k] )
        call c_f_pointer( jmodel%data, jdata, [k] )
        geometric = sqrt(idata(k)*jdata(k))
      end function geometric
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end function cross_pair

!===================================================================================================
!                                      B O N D     M O D E L S
!===================================================================================================

  function EmDee_bond_harmonic( k, r0 ) result(model) bind(C,name="EmDee_bond_harmonic")
    real(rb), value :: k, r0
    type(EmDee_Model) :: model
    model = EmDee_Model( mHARMOMIC, set_data( [k, r0] ), r0, -k, 0.5_rb*k )
  end function EmDee_bond_harmonic

!---------------------------------------------------------------------------------------------------

  function EmDee_bond_morse( D, alpha, r0 ) result(model) bind(C,name="EmDee_bond_morse")
    real(rb), value :: D, alpha, r0
    type(EmDee_Model) :: model
    model = EmDee_Model( mMORSE, set_data( [D, alpha, r0] ), r0, -alpha, D, -2.0_rb*D*alpha )
  end function EmDee_bond_morse

!===================================================================================================
!                                    A N G L E     M O D E L S
!===================================================================================================

  function EmDee_angle_harmonic( k, theta0 ) result(model) bind(C,name="EmDee_angle_harmonic")
    real(rb), value :: k, theta0
    type(EmDee_Model) :: model
    model = EmDee_Model( mHARMOMIC, set_data( [k, theta0] ), Deg2Rad*theta0, -k, 0.5_rb*k )
  end function EmDee_angle_harmonic

!===================================================================================================
!                                 D I H E D R A L     M O D E L S
!===================================================================================================

  function EmDee_dihedral_harmonic( k, phi0 ) result(model) bind(C,name="EmDee_dihedral_harmonic")
    real(rb), value :: k, phi0
    type(EmDee_Model) :: model
    model = EmDee_Model( mHARMOMIC, set_data( [k, phi0] ), 0.0_rb, Deg2Rad*phi0, -k, 0.5_rb*k )
  end function EmDee_dihedral_harmonic

!===================================================================================================
!                             A U X I L I A R Y     P R O C E D U R E
!===================================================================================================

  function set_data( values ) result( data )
    real(rb),     intent(in) :: values(:)
    type(c_ptr)              :: data
    real(rb), pointer :: ptr(:)
    allocate( ptr(size(values)) )
    ptr = values
    data = c_loc(ptr)
  end function set_data

!===================================================================================================

end module models
