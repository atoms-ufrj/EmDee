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

! INPUT:  invR2, Qi, Qj, model
! OUTPUT: Eij, Wij

block
  real(rb) :: rFc, invR, QiQj, QiQjbyR
  invR = sqrt(invR2)
  QiQj = Qi*Qj
  QiQjbyR = QiQj*invR
  rFc = QiQj*model%fshift_coul/invR
  Eij = QiQjbyR + QiQj*model%eshift_coul + rFc
  Wij = QiQjbyR - rFc
end block
