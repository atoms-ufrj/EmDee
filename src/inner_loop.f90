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

  if (r2 < Rc2) then
    invR2 = me%invL2/r2
    jtype = me%atomType(j)
    Qj = me%charge(j)
    if (compute) then
      associate( model => partner(jtype)%model )
        select type ( model )
          include "compute_pair.f90"
        end select
        if (icharged.and.me%charged(j)) then
          if (model%noInvR) invR = sqrt(invR2)
          QiQj = Qi*Qj
          QiQjbyR = QiQj*invR
          rFc = QiQj*me%fshift/invR
          Eij = Eij + QiQjbyR + QiQj*me%eshift + rFc
          Wij = Wij + QiQjbyR - rFc
        end if
      end associate
      Potential = Potential + Eij
      Virial = Virial + Wij
      Fij = Wij*invR2*Rij
      if (multilayer(jtype)) then
        associate( pair => me%pair(itype,jtype,:) )
          do layer = 1, me%nlayers
            select type ( model => pair(layer)%model )
              include "compute_pair.f90"
            end select
            Elayer(layer) = Elayer(layer) + Eij
            Wlayer(layer) = Wlayer(layer) + Wij
          end do
        end associate
      end if
    else
      associate( model => partner(jtype)%model )
        select type ( model )
          include "virial_compute_pair.f90"
        end select
        if (icharged.and.me%charged(j)) then
          if (model%noInvR_virial) invR = sqrt(invR2)
          Wij = Wij + Qi*Qj*(invR - me%fshift/invR)
        end if
      end associate
      Virial = Virial + Wij
      Fij = Wij*invR2*Rij
    end if
    Fi = Fi + Fij
    F(:,j) = F(:,j) - Fij
  end if

