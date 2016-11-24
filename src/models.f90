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

use pair_lj_cut_module
use pair_lj_cut_coul_cut_module
use pair_lj_sf_module
use pair_lj_sf_coul_sf_module

implicit none

integer, parameter :: mCOULOMB = 10000

integer, parameter :: mNONE     = 0, &  ! No model
                      mLJ       = 1, &  ! Lennard-Jones
                      mLJSF     = 2, &  ! Lennard-Jones, Shifted-Force
                      mSOFTCORE = 3, &  ! Softcore Model, Beutler et al. (1994)
                      mHARMOMIC = 4, &  ! Harmonic
                      mMORSE    = 5     ! Morse

real(rb), parameter, private :: Deg2Rad = 3.14159265358979324_rb/180_rb

type, extends(cModel) :: tModel
  real(rb) :: p1 = zero
  real(rb) :: p2 = zero
  real(rb) :: p3 = zero
  real(rb) :: p4 = zero
  real(rb) :: factor = zero
  contains
    procedure :: compute => tModel_compute
end type tModel

contains

!---------------------------------------------------------------------------------------------------

  type(c_ptr) function EmDee_model_none() bind(C,name="EmDee_model_none")
    type(tModel), pointer :: model
    allocate( model )
    model%id = mNONE
    EmDee_model_none = c_loc(model)
  end function EmDee_model_none

!===================================================================================================
!                                      P A I R     M O D E L S
!===================================================================================================

!  function EmDee_pair_lj( epsilon, sigma ) bind(C,name="EmDee_pair_lj")
!    real(rb), value :: epsilon, sigma
!    type(c_ptr)     :: EmDee_pair_lj

!    type(cModelPtr), pointer :: container

!    allocate( container )
!    allocate( pair_lj::container%model )

!    select type (model => container%model)
!      type is (pair_lj)
!        model%id = mLJ
!        model%data => set_data( [epsilon, sigma] )
!        model%p1 = sigma*sigma
!        model%p2 = 4.0_rb*epsilon

!        model%epsilon = epsilon
!        model%sigma = sigma
!        model%eps4 = 4.0_rb*epsilon
!        model%sigsq = sigma*sigma
!    end select
!    EmDee_pair_lj = c_loc(container)

!  end function EmDee_pair_lj

  subroutine tModel_compute( model, Eij, Wij, invR2, Qi, Qj )
    class(tModel), intent(in)  :: model
    real(rb),       intent(out) :: Eij, Wij
    real(rb),       intent(in)  :: invR2, Qi, Qj
    Eij = zero
    Wij = zero
  end subroutine tModel_compute

!---------------------------------------------------------------------------------------------------

!  function EmDee_pair_lj_sf_old( epsilon, sigma, cutoff ) bind(C,name="EmDee_pair_lj_sf_old")
!    real(rb), value :: epsilon, sigma, cutoff
!    type(c_ptr)     :: EmDee_pair_lj_sf_old

!    type(tModel), pointer :: model
!    real(rb) :: sr6, sr12, eps4, Ec, Fc

!    sr6 = (sigma/cutoff)**6
!    sr12 = sr6*sr6
!    eps4 = 4.0_rb*epsilon
!    Ec = eps4*(sr12 - sr6)
!print*, "Evdw = ", Ec
!    Fc = 6.0_rb*(eps4*sr12 + Ec)/cutoff
!print*, "Wvdw = ", Fc*cutoff
!    Ec = -(Ec + Fc*cutoff)
!    allocate( model )
!    model = tModel( mLJSF, set_data( [epsilon, sigma, cutoff] ), sigma**2, eps4, Fc, Ec )
!    EmDee_pair_lj_sf_old = c_loc(model)

!print*, Ec, Fc

!  end function EmDee_pair_lj_sf_old

!---------------------------------------------------------------------------------------------------

!  function EmDee_pair_softcore( epsilon, sigma, lambda ) bind(C,name="EmDee_pair_softcore")
!    real(rb), value :: epsilon, sigma, lambda
!    type(c_ptr)     :: EmDee_pair_softcore

!    type(tModel), pointer :: model

!    allocate( model )
!    model = tModel( mSOFTCORE, set_data( [epsilon, sigma, lambda] ), sigma**2 )
!    EmDee_pair_softcore = c_loc(model)

!  end function EmDee_pair_softcore

!===================================================================================================
!                                     M I X I N G     R U L E S
!===================================================================================================

!  function cross_pair( i, j ) result( ij )
!    class(cPairModel), pointer, intent(in) :: i, j
!    class(cPairModel), pointer             :: ij

!    type(c_ptr)              :: ijmodel
!    type(pairModelContainer), pointer :: container

!    if (associated(i).and.associated(j)) then

!      if (match(mLJ,mLJ)) then
!        ijmodel = EmDee_pair_lj( epsilon = geometric(1), &
!                                 sigma = arithmetic(2)   )

!      else if (match(mLJSF,mLJSF)) then
!        ijmodel = EmDee_pair_lj_sf( epsilon = geometric(1), &
!                                    sigma = arithmetic(2),  &
!                                    cutoff = arithmetic(3)  )

!      else if (match(mSOFTCORE,mSOFTCORE)) then
!        ijmodel = EmDee_pair_lj( epsilon = geometric(1), &
!                                 sigma = arithmetic(2)   )

!      else if (match(mSOFTCORE,mLJ)) then
!        ijmodel = EmDee_pair_softcore( epsilon = geometric(1),    &
!                                       sigma = arithmetic(2),     &
!                                       lambda = from(mSOFTCORE,3) )

!      else if ((i%id == mNONE).or.(j%id == mNONE)) then
!        ijmodel = EmDee_model_none()

!      else
!        ijmodel = c_null_ptr

!      end if
!    else
!      ijmodel = c_null_ptr
!    end if

!    if (c_associated(ijmodel)) then
!      call c_f_pointer( ijmodel, container )
!!      ij => container%model
!!      ij%external = .false.
!    end if

!    contains
!      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      logical function match( a, b )
!        integer, intent(in) :: a, b
!        match = ((i%id == a).and.(j%id == b)).or.((i%id == b).and.(j%id == a))
!      end function match
!      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      real(rb) function arithmetic( k )
!        integer, intent(in) :: k
!        arithmetic = 0.5_rb*(i%data(k) + j%data(k))
!      end function arithmetic
!      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      real(rb) function geometric( k )
!        integer, intent(in) :: k
!        geometric = sqrt(i%data(k)*j%data(k))
!      end function geometric
!      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      real(rb) function from( id, k )
!        integer, intent(in) :: id, k
!        if (id == i%id) then
!          from = i%data(k)
!        else if (id == j%id) then
!          from = j%data(k)
!        else
!          stop "ERROR defining mixing rule"
!        end if
!      end function from
!      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  end function cross_pair

!===================================================================================================
!                                      B O N D     M O D E L S
!===================================================================================================

  type(c_ptr) function EmDee_bond_harmonic( k, r0 ) bind(C,name="EmDee_bond_harmonic")
    real(rb), value :: k, r0
    type(tModel), pointer :: model
    allocate( model )
    model = tModel( mHARMOMIC, set_data( [k, r0] ), r0, -k, 0.5_rb*k )
    EmDee_bond_harmonic = c_loc(model)
  end function EmDee_bond_harmonic

!---------------------------------------------------------------------------------------------------

!  type(c_ptr) function EmDee_bond_morse( D, alpha, r0 ) bind(C,name="EmDee_bond_morse")
!    real(rb), value :: D, alpha, r0
!    type(tModel), pointer :: model
!    allocate( model )
!    model = tModel( mMORSE, set_data( [D, alpha, r0] ), r0, -alpha, D, -2.0_rb*D*alpha )
!    EmDee_bond_morse = c_loc(model)
!  end function EmDee_bond_morse

!===================================================================================================
!                                    A N G L E     M O D E L S
!===================================================================================================

!  type(c_ptr) function EmDee_angle_harmonic( k, theta0 ) bind(C,name="EmDee_angle_harmonic")
!    real(rb), value :: k, theta0
!    type(tModel), pointer :: model
!    allocate( model )
!    model = tModel( mHARMOMIC, set_data( [k, theta0] ), Deg2Rad*theta0, -k, 0.5_rb*k )
!    EmDee_angle_harmonic = c_loc(model)
!  end function EmDee_angle_harmonic

!===================================================================================================
!                                 D I H E D R A L     M O D E L S
!===================================================================================================

!  type(c_ptr) function EmDee_dihedral_harmonic( k, phi0 ) bind(C,name="EmDee_dihedral_harmonic")
!    real(rb), value :: k, phi0
!    type(tModel), pointer :: model
!    allocate( model )
!    model = tModel( mHARMOMIC, set_data( [k, phi0] ), 0.0_rb, Deg2Rad*phi0, -k, 0.5_rb*k )
!    EmDee_dihedral_harmonic = c_loc(model)
!  end function EmDee_dihedral_harmonic

!===================================================================================================
!                             A U X I L I A R Y     P R O C E D U R E
!===================================================================================================

  function set_data( values ) result( data )
    real(rb), intent(in) :: values(:)
    real(rb), pointer    :: data(:)
    allocate( data(size(values)) )
    data = values
  end function set_data

!===================================================================================================

end module models

