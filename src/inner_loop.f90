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

block

  integer  :: jtype
  real(rb) :: invR2, invR, Wij, QiQj, WCij, rFc
  real(rb) :: Rij(3), Fij(3)
  logical  :: ijcharged, noInvR
# if defined(compute)
  integer  :: l, layer
  real(rb) :: Eij, ECij
# endif

  Rij = Rvec(:,m)
  r2 = sum(Rij*Rij)
  if (r2 < Rc2) then
    invR2 = me%invL2/r2
    invR = sqrt(invR2)
    jtype = me%atomType(j)
    ijcharged = icharged.and.me%charged(j)
    noInvR = .false.
    associate( pair => partner(jtype) )
      associate ( model => pair%model )
        select type ( model )
#         if defined(compute)
            include "compute_pair.f90"
#         else
            include "virial_compute_pair.f90"
#         endif
        end select
        if (model%shifted_force) then
          if (noInvR) then
            invR = sqrt(invR2)
            noInvR = .false.
          end if
          rFc = model%fshift/invR
          Wij = Wij - rFc
#         if defined(compute)
            Eij = Eij + model%eshift + rFc
#         endif
        end if
      end associate
#     if defined(compute)
        Epair = Epair + Eij
#     endif
      Wpair = Wpair + Wij
      if (ijcharged.and.pair%coulomb) then
        QiQj = pair%kCoul*Qi*me%charge(j)
        select type ( model => me%coul(me%layer)%model )
#         if defined(compute)
            include "compute_coul.f90"
#         else
            include "virial_compute_coul.f90"
#         endif
        end select
#       if defined(compute)
          Ecoul = Ecoul + ECij
          Eij = Eij + ECij
#       endif
        Wcoul = Wcoul + WCij
        Wij = Wij + WCij
      end if
    end associate
    Fij = Wij*invR2*Rij
    Fi = Fi + Fij
    F(:,j) = F(:,j) - Fij

#   if defined(compute)
      if (multilayer(jtype)) then
        Elayer(me%layer) = Elayer(me%layer) + Eij
        do l = 1, me%nlayers-1
          layer = me%other_layer(l)
          associate( pair => me%pair(itype,jtype,layer) )
            select type ( model => pair%model )
              include "energy_compute_pair.f90"
            end select
            if (ijcharged.and.pair%coulomb) then
              QiQj = pair%kCoul*Qi*me%charge(j)
              select type ( model => me%coul(me%layer)%model )
                include "energy_compute_coul.f90"
              end select
              Eij = Eij + ECij
            end if
            Elayer(layer) = Elayer(layer) + Eij
          end associate
        end do
      end if
#   endif
  end if

end block
