!---------------------------------------------------------------------------------------------------
subroutine compute_pair
  select case (ij%model)
    case (mLJ)
      call lennard_jones_compute( Eij, Wij, invR2*ij%p1, ij%p2 )
    case (mSFLJ)
      call lennard_jones_sf_compute( Eij, Wij, invR2*ij%p1, ij%p2, ij%p3/sqrt(invR2), ij%p4 )
  end select
end subroutine compute_pair
!---------------------------------------------------------------------------------------------------
pure subroutine lennard_jones_compute( E, W, sr2, eps4 )
  real(rb), intent(out) :: E, W
  real(rb), intent(in)  :: sr2, eps4
  real(rb) :: sr6, sr12
  sr6 = sr2*sr2*sr2
  sr12 = sr6*sr6
  E = eps4*(sr12 - sr6)
  W = 6.0_rb*(eps4*sr12 + E)
end subroutine lennard_jones_compute
!---------------------------------------------------------------------------------------------------
pure subroutine lennard_jones_sf_compute( E, W, sr2, eps4, rFc, shift )
  real(rb), intent(out) :: E, W
  real(rb), intent(in)  :: sr2, eps4, rFc, shift
  real(rb) :: sr6, sr12
  sr6 = sr2*sr2*sr2
  sr12 = sr6*sr6
  E = eps4*(sr12 - sr6)
  W = 6.0_rb*(eps4*sr12 + E) - rFc
  E = E + rFc + shift
end subroutine lennard_jones_sf_compute
!---------------------------------------------------------------------------------------------------

