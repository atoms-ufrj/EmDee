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
  real(rb) :: sr2, sr6, sr12, rFc, invR, QiQj, QiQjbyR
  sr2 = model%sigSq*invR2
  sr6 = sr2*sr2*sr2
  sr12 = sr6*sr6
  Eij = model%eps4*(sr12 - sr6)
  invR = sqrt(invR2)
  QiQj = Qi*Qj
  QiQjbyR = QiQj*invR
  rFc = QiQjbyR*model%fshift_coul
  Wij = 6.0_rb*(model%eps4*sr12 + Eij) + QiQjbyR - rFc
  Eij = Eij + QiQjbyR + rFc + QiQj*model%eshift_coul
end block