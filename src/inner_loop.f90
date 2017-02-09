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
  ijcharged = icharged.and.me%charged(j)
  noInvR = .true.
  if (compute) then
    associate( pair => partner(jtype) )
      select type ( model => pair%model )
        include "compute_pair.f90"
      end select
      if (ijcharged.and.pair%coulomb) then
        QiQj = pair%kCoul*Qi*me%charge(j)
        select type ( model => me%coul(me%layer)%model )
          include "compute_coul.f90"
        end select
      else
        ECij = zero
        WCij = zero
      end if
    end associate
    Wij = Wij + WCij
    Virial = Virial + Wij
    Fij = Wij*invR2*Rij
    Epair = Epair + Eij
    Ecoul = Ecoul + ECij
    if (multilayer(jtype)) then
      Elayer(me%layer) = Elayer(me%layer) + Eij
      Wlayer(me%layer) = Wlayer(me%layer) + Wij
      associate( pair => me%pair(itype,jtype,:) )
        do l = 1, me%nlayers-1
          layer = me%other_layer(l)
          select type ( model => pair(layer)%model )
            include "compute_pair.f90"
          end select
          Elayer(layer) = Elayer(layer) + Eij
          Wlayer(layer) = Wlayer(layer) + Wij
        end do
      end associate
    end if
  else
    associate( pair => partner(jtype) )
      select type ( model => pair%model )
        include "virial_compute_pair.f90"
      end select
      if (ijcharged.and.pair%coulomb) then
        QiQj = pair%kCoul*Qi*me%charge(j)
        select type ( model => me%coul(me%layer)%model )
          include "virial_compute_coul.f90"
        end select
      end if
    end associate
    Virial = Virial + Wij
    Fij = Wij*invR2*Rij
  end if
  Fi = Fi + Fij
  F(:,j) = F(:,j) - Fij
end if

