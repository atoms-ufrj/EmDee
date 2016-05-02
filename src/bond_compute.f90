!---------------------------------------------------------------------------------------------------
subroutine compute_bond( ij, r )
  type(tModel), intent(in) :: ij
  real(rb),     intent(in) :: r
  select case (ij%model)
    case (HARM)
      call harmonic_compute( E, F, r - ij%p1, ij%p2, ij%p3 )
  end select
end subroutine compute_bond
!---------------------------------------------------------------------------------------------------
pure subroutine harmonic_compute( E, F, x_minus_x0, k, k_2 )
  real(rb), intent(out) :: E, F
  real(rb), intent(in)  :: x_minus_x0, k, k_2
  E = k_2*x_minus_x0*x_minus_x0
  F = k*x_minus_x0
end subroutine harmonic_compute
!---------------------------------------------------------------------------------------------------

