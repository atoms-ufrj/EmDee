subroutine compute_pair
  Rij = Ri - Rs(:,j)
  Rij = Rij - nint(Rij)
  r2 = sum(Rij*Rij)
  if (r2 < Rc2) then
    invR2 = invL2/r2
    ij => pairType(itype,type(j))
    select case (ij%model)
      case (LJ)
        call lennard_jones( Eij, Wij, invR2*ij%p1, ij%p2 )
      case (SF_LJ)
        call lennard_jones_sf( Eij, Wij, invR2*ij%p1, ij%p2, sqrt(r2)*ij%p3, ij%p4 )
    end select
    Epot = Epot + Eij
    Virial = Virial + Wij
    Fij = Wij*invR2*Rij
    Fi = Fi + Fij
    Fs(:,j) = Fs(:,j) - Fij
  end if
end subroutine compute_pair

subroutine lennard_jones( E, W, sr2, eps4 )
  real(rb), intent(out) :: E, W
  real(rb), intent(in)  :: sr2, eps4
  real(rb) :: sr6, sr12
  sr6 = sr2*sr2*sr2
  sr12 = sr6*sr6
  E = eps4*(sr12 - sr6)
  W = 6.0_rb*(eps4*sr12 + E)
end subroutine lennard_jones
!---------------------------------------------------------------------------------------------------
subroutine lennard_jones_sf( E, W, sr2, eps4, rFc, shift )
  real(rb), intent(out) :: E, W
  real(rb), intent(in)  :: sr2, eps4, rFc, shift
  real(rb) :: sr6, sr12
  sr6 = sr2*sr2*sr2
  sr12 = sr6*sr6
  E = eps4*(sr12 - sr6)
  W = 6.0_rb*(eps4*sr12 + E) - rFc
  E = E + rFc + shift
end subroutine lennard_jones_sf

