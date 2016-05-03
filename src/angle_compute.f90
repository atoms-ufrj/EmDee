!---------------------------------------------------------------------------------------------------
subroutine compute_angle( ijk, theta )
  type(tModel), intent(in) :: ijk
  real(rb),     intent(in) :: theta
  select case (ijk%model)
    case (mHARMOMIC)
      call harmonic_compute( E, F, theta - ijk%p1, ijk%p2, ijk%p3 )
  end select
end subroutine compute_angle
!---------------------------------------------------------------------------------------------------
pure subroutine harmonic_compute( E, F, x_minus_x0, k, k_2 )
  real(rb), intent(out) :: E, F
  real(rb), intent(in)  :: x_minus_x0, k, k_2
  E = k_2*x_minus_x0*x_minus_x0
  F = k*x_minus_x0
end subroutine harmonic_compute
!---------------------------------------------------------------------------------------------------

