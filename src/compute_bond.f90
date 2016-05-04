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

!---------------------------------------------------------------------------------------------------

subroutine compute_bond
  select case (model%id)
    case (mHARMOMIC)
      call harmonic( E, mdEdr, d - model%p1, model%p2, model%p3 )
    case (mMORSE)
      call morse( E, mdEdr, exp(model%p2*(d - model%p1)), model%p2, model%p3 )
  end select
end subroutine compute_bond

!---------------------------------------------------------------------------------------------------

pure subroutine harmonic( E, F, x_minus_x0, minus_k, k_2 )
  real(rb), intent(out) :: E, F
  real(rb), intent(in)  :: x_minus_x0, minus_k, k_2
  E = k_2*x_minus_x0*x_minus_x0
  F = minus_k*x_minus_x0
end subroutine harmonic

!---------------------------------------------------------------------------------------------------

pure subroutine morse( E, F, expAdeltaR, D, m2Dalpha )
  real(rb), intent(out) :: E, F
  real(rb), intent(in)  :: expAdeltaR, D, m2Dalpha
  real(rb) :: x
  x = 1.0_rb - expAdeltaR
  E = D*x*x
  F = m2Dalpha*x*expAdeltaR
end subroutine morse

!---------------------------------------------------------------------------------------------------

