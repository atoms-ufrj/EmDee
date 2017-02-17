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

Rij = Rvec(:,m)
r2 = sum(Rij*Rij)
if (r2 < Rc2) then
  invR2 = me%invL2/r2
  jtype = me%atomType(j)
  ijcharged = icharged.and.me%charged(j)
  noInvR = .true.
  associate( pair => partner(jtype) )
    select type ( model => pair%model )
      include "virial_compute_pair.f90"
    end select
    Wpair = Wpair + Wij
    if (ijcharged.and.pair%coulomb) then
      QiQj = pair%kCoul*Qi*me%charge(j)
      select type ( model => me%coul(me%layer)%model )
        include "virial_compute_coul.f90"
      end select
      Wcoul = Wcoul + WCij
      Wij = Wij + WCij
    end if
  end associate
  Fij = Wij*invR2*Rij
  Fi = Fi + Fij
  F(:,j) = F(:,j) - Fij
end if
