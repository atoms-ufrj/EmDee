!---------------------------------------------------------------------------------------------------

subroutine compute_bond( ij, r )
  type(tModel), intent(in) :: ij
  real(rb),     intent(in) :: r
  select case (ij%model)
    case (mHARMOMIC)
      call harmonic_compute( E, F, r - ij%p1, ij%p2, ij%p3 )
    case (mMORSE)
      call morse_compute( E, F, exp(ij%p2*(r - ij%p1)), ij%p2, ij%p3 )
  end select
end subroutine compute_bond

!---------------------------------------------------------------------------------------------------

pure subroutine harmonic_compute( E, F, x_minus_x0, minus_k, k_2 )
  real(rb), intent(out) :: E, F
  real(rb), intent(in)  :: x_minus_x0, minus_k, k_2
  E = k_2*x_minus_x0*x_minus_x0
  F = minus_k*x_minus_x0
end subroutine harmonic_compute

!---------------------------------------------------------------------------------------------------

pure subroutine morse_compute( E, F, expAdeltaR, D, m2Dalpha )
  real(rb), intent(out) :: E, F
  real(rb), intent(in)  :: expAdeltaR, D, m2Dalpha
  real(rb) :: x
  x = 1.0_rb - expAdeltaR
  E = D*x*x
  F = m2Dalpha*x*expAdeltaR
end subroutine morse_compute

!---------------------------------------------------------------------------------------------------

