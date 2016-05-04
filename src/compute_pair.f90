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

subroutine compute_pair()
  select case (model%id)
    case (mLJ)
      call lj( Eij, Wij, invR2*model%p1, model%p2 )
    case (mLJSF)
      call lj_sf( Eij, Wij, invR2*model%p1, model%p2, model%p3/sqrt(invR2), model%p4 )
    case (mLJCOULSF)
      call lj_coul_sf( Eij, Wij, invR2*model%p1, model%p2, model%p3/sqrt(invR2), model%p4 )
  end select
end subroutine compute_pair

!---------------------------------------------------------------------------------------------------

pure subroutine lj( E, W, sr2, eps4 )
  real(rb), intent(out) :: E, W
  real(rb), intent(in)  :: sr2, eps4
  real(rb) :: sr6, sr12
  sr6 = sr2*sr2*sr2
  sr12 = sr6*sr6
  E = eps4*(sr12 - sr6)
  W = 6.0_rb*(eps4*sr12 + E)
end subroutine lj

!---------------------------------------------------------------------------------------------------

pure subroutine lj_sf( E, W, sr2, eps4, rFc, shift )
  real(rb), intent(out) :: E, W
  real(rb), intent(in)  :: sr2, eps4, rFc, shift
  real(rb) :: sr6, sr12
  sr6 = sr2*sr2*sr2
  sr12 = sr6*sr6
  E = eps4*(sr12 - sr6)
  W = 6.0_rb*(eps4*sr12 + E) - rFc
  E = E + rFc + shift
end subroutine lj_sf

!---------------------------------------------------------------------------------------------------

pure subroutine lj_coul_sf( E, W, sr2, eps4, rFc, shift )
  real(rb), intent(out) :: E, W
  real(rb), intent(in)  :: sr2, eps4, rFc, shift
  real(rb) :: sr6, sr12
  sr6 = sr2*sr2*sr2
  sr12 = sr6*sr6
  E = eps4*(sr12 - sr6)
  W = 6.0_rb*(eps4*sr12 + E) - rFc
  E = E + rFc + shift
end subroutine lj_coul_sf

!---------------------------------------------------------------------------------------------------

