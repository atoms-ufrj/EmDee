!---------------------------------------------------------------------------------------------------
type(tPairType) function lennard_jones( sigma, epsilon ) bind(C)
  real(rb), value :: sigma, epsilon

  lennard_jones = tPairType( LJ, sigma*sigma, 4.0_rb*epsilon, 0.0_rb, 0.0_rb )

end function lennard_jones
!---------------------------------------------------------------------------------------------------
type(tPairType) function shifted_force_lennard_jones( sigma, epsilon, rc ) bind(C)
  real(rb), value :: sigma, epsilon, rc

  real(rb) :: sr6, sr12, eps4, Ec, Fc
  sr6 = (sigma/rc)**6
  sr12 = sr6*sr6
  eps4 = 4.0_rb*epsilon
  Ec = eps4*(sr12 - sr6)
  Fc = 6.0_rb*(eps4*sr12 + Ec)/rc
  shifted_force_lennard_jones = tPairType( SF_LJ, sigma**2, eps4, Fc, -(Ec + Fc*rc) )

end function shifted_force_lennard_jones
!---------------------------------------------------------------------------------------------------

