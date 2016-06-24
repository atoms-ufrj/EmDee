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
                          mLJSF     = 2, &  ! Lennard-Jones (Shifted-Force)
                          mLJCOULSF = 3, &  ! Lennard-Jones/Coulomb (Shifted-Force)
                          mHARMOMIC = 4, &  ! Harmonic
                          mMORSE    = 5     ! Morse

type, bind(C) :: EmDee_Model
  type(c_ptr) :: data
  type(c_ptr) :: params
  integer(ib) :: external = 1
end type EmDee_Model

type md_data
  real(rb), allocatable :: value(:)
end type md_data

type md_params
  integer(ib) :: id = mNONE
  real(rb)    :: p1 = 0.0_rb
  real(rb)    :: p2 = 0.0_rb
  real(rb)    :: p3 = 0.0_rb
  real(rb)    :: p4 = 0.0_rb
end type md_params

type param_ptr
  type(md_params), pointer :: params => null()
end type param_ptr

type model_ptr
  type(EmDee_Model), pointer :: model => null()
end type model_ptr

private :: set_data, set_params

contains

!===================================================================================================
!                                      P A I R     M O D E L S
!===================================================================================================

  function EmDee_pair_lj( sigma, epsilon ) result(model) bind(C,name="EmDee_pair_lj")
    real(rb), value   :: sigma, epsilon
    type(EmDee_Model) :: model

    model%data = set_data( [sigma, epsilon] )
    model%params = set_params( mLJ, sigma*sigma, 4.0_rb*epsilon )

  end function EmDee_pair_lj

!---------------------------------------------------------------------------------------------------

  function EmDee_pair_lj_sf( sigma, epsilon, cutoff ) result( model ) &
                                                        bind(C,name="EmDee_pair_lj_sf")
    real(rb), value :: sigma, epsilon, cutoff
    type(EmDee_Model) :: model

    real(rb) :: sr6, sr12, eps4, Ec, Fc

    sr6 = (sigma/cutoff)**6
    sr12 = sr6*sr6
    eps4 = 4.0_rb*epsilon
    Ec = eps4*(sr12 - sr6)
    Fc = 6.0_rb*(eps4*sr12 + Ec)/cutoff

    model%data = set_data( [sigma, epsilon, cutoff] )
    model%params = set_params( mLJSF, sigma**2, eps4, Fc, -(Ec + Fc*cutoff) )

  end function EmDee_pair_lj_sf

!===================================================================================================
!                                     M I X I N G     R U L E S
!===================================================================================================

  function cross_pair( imodel, jmodel ) result( ij )
    type(EmDee_Model), pointer, intent(in) :: imodel, jmodel
    type(EmDee_Model), pointer             :: ij

    type(md_data),   pointer :: idata, jdata
    type(md_params), pointer :: iparams, jparams

    if (associated(imodel).and.associated(jmodel)) then
      call c_f_pointer( imodel%data, idata )
      call c_f_pointer( jmodel%data, jdata )
      call c_f_pointer( imodel%params, iparams )
      call c_f_pointer( jmodel%params, jparams )

      allocate( ij )
      if (match(mLJ,mLJ)) then
        ij = EmDee_pair_lj( arithmetic(1), geometric(2) )
      else if (match(mLJSF,mLJSF).or.match(mLJSF,mLJCOULSF)) then
        ij = EmDee_pair_lj_sf( arithmetic(1), geometric(2), arithmetic(3) )
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
        match = ((iparams%id == a).and.(jparams%id == b)) .or. &
                ((iparams%id == b).and.(jparams%id == a))
      end function match
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
  end function cross_pair

!===================================================================================================
!                                      B O N D     M O D E L S
!===================================================================================================

  function EmDee_bond_harmonic( k, r0 ) result( model ) bind(C,name="EmDee_bond_harmonic")
    real(rb), value :: k, r0
    type(EmDee_Model) :: model

    model%data = set_data( [k, r0] )
    model%params = set_params( mHARMOMIC, r0, -k, 0.5_rb*k )

  end function EmDee_bond_harmonic

!---------------------------------------------------------------------------------------------------

  function EmDee_bond_morse( D, alpha, r0 ) result( model ) bind(C,name="EmDee_bond_morse")
    real(rb), value :: D, alpha, r0
    type(EmDee_Model) :: model

    model%data = set_data( [D, alpha, r0] )
    model%params = set_params( mMORSE, r0, -alpha, D, -2.0_rb*D*alpha )

  end function EmDee_bond_morse

!===================================================================================================
!                                    A N G L E     M O D E L S
!===================================================================================================

  function EmDee_angle_harmonic( k, theta0 ) result( model ) bind(C,name="EmDee_angle_harmonic")
    real(rb), value :: k, theta0
    type(EmDee_Model) :: model

    model%data = set_data( [k, theta0] )
    model%params = set_params( mHARMOMIC, theta0, -k, 0.5_rb*k )

  end function EmDee_angle_harmonic

!===================================================================================================
!                                 D I H E D R A L     M O D E L S
!===================================================================================================

  function EmDee_dihedral_harmonic( k, phi0 ) result( model ) bind(C,name="EmDee_dihedral_harmonic")
    real(rb), value :: k, phi0
    type(EmDee_Model) :: model

    model%data = set_data( [k, phi0] )
    model%params = set_params( mHARMOMIC, 0.0_rb, phi0, -k, 0.5_rb*k )

  end function EmDee_dihedral_harmonic

!===================================================================================================
!                             A U X I L I A R Y     P R O C E D U R E S
!===================================================================================================

  function set_data( values ) result( data )
    real(rb),     intent(in) :: values(:)
    type(c_ptr)              :: data

    type(md_data), pointer :: ptr

    allocate( ptr )
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
