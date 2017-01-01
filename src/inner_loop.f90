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
    select type ( model => partner(jtype)%model )
      include "compute_pair.f90"
    end select
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
    select type ( model => partner(jtype)%model )
      include "compute_pair_virial.f90"
    end select
    Fij = Wij*invR2*Rij
  end if
  Fi = Fi + Fij
  F(:,j) = F(:,j) - Fij
end if
