module model_setup

use, intrinsic :: iso_c_binding

implicit none

integer, parameter, private :: ib = c_int, &
                               rb = c_double

type, bind(C) :: tModel
  integer(ib) :: model
  real(rb)    :: p1, p2, p3, p4
end type tModel

integer(ib), parameter :: mNONE  = 0, &
                          mLJ    = 1, &
                          mSFLJ = 2, &
                          mHARMOMIC = 3, &
                          mMORSE    = 4

contains

!---------------------------------------------------------------------------------------------------

type(tModel) function lennard_jones( sigma, epsilon ) bind(C)
  real(rb), value :: sigma, epsilon
  lennard_jones = tModel( mLJ, sigma*sigma, 4.0_rb*epsilon, 0.0_rb, 0.0_rb )
end function lennard_jones

!---------------------------------------------------------------------------------------------------

type(tModel) function shifted_force_lennard_jones( sigma, epsilon, rc ) bind(C)
  real(rb), value :: sigma, epsilon, rc
  real(rb) :: sr6, sr12, eps4, Ec, Fc
  sr6 = (sigma/rc)**6
  sr12 = sr6*sr6
  eps4 = 4.0_rb*epsilon
  Ec = eps4*(sr12 - sr6)
  Fc = 6.0_rb*(eps4*sr12 + Ec)/rc
  shifted_force_lennard_jones = tModel( mSFLJ, sigma**2, eps4, Fc, -(Ec + Fc*rc) )
end function shifted_force_lennard_jones

!---------------------------------------------------------------------------------------------------

type(tModel) function harmonic( k, x0 ) bind(C)
  real(rb), value :: k, x0
  harmonic = tModel( mHARMOMIC, x0, -k, 0.5_rb*k, 0.0_rb )
end function harmonic

!---------------------------------------------------------------------------------------------------

type(tModel) function morse( D, alpha, r0 ) bind(C)
  real(rb), value :: D, alpha, r0
  morse = tModel( mMORSE, r0, -alpha, D, -2.0_rb*D*alpha )
end function morse

!---------------------------------------------------------------------------------------------------

end module model_setup
