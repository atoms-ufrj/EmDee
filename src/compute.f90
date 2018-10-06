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
  integer  :: i, j, k, m, itype, jtype, firstAtom, lastAtom
  real(rb) :: r2, invR2, invR, Wij, Wsum, Qi, QiQj, rFc, Eij
  real(rb) :: Rij(3), Ri(3), Fi(3), Fij(3)
  logical  :: icharged, ijcharged
  real(rb), allocatable :: Rvec(:,:)

  associate ( neighbor => me%neighbor(thread) )
    firstAtom = me%cellAtom%first(me%threadCell%first(thread))
    lastAtom = me%cellAtom%last(me%threadCell%last(thread))
    do k = firstAtom, lastAtom
      i = me%cellAtom%item(k)
      itype = me%atomType(i)
      Qi = me%charge(i)
      icharged = me%charged(i)
      Ri = Rs(:,i)
      Fi = zero
      associate ( partner => me%pair(:,itype,me%layer), &
                  jlist => neighbor%item(neighbor%first(i):upper(i)) )
        Rvec = Rs(:,jlist)
        forall (m=1:size(jlist)) Rvec(:,m) = pbc(Ri - Rvec(:,m))
        do m = 1, size(jlist)
          j = jlist(m)
          Rij = Rvec(:,m)
          r2 = sum(Rij*Rij)
          if (r2 < Rc2) then
            invR2 = invL2/r2
            invR = sqrt(invR2)
            jtype = me%atomType(j)
            ijcharged = icharged.and.me%charged(j)
            associate( pair => partner(jtype) )
              associate ( model => pair%model )
                select type ( model )
#                 if defined(compute)
#                   include "compute_pair.f90"
#                 else
#                   include "virial_compute_pair.f90"
#                 endif
                end select
#               include "apply_modifier.f90"
              end associate
#             if defined(compute)
                Epair = Epair + Eij
#             endif
              Wpair = Wpair + Wij
              Wsum = Wij
              if (ijcharged.and.pair%coulomb) then
                QiQj = pair%kCoul*Qi*me%charge(j)
                select type ( model => me%coul(me%layer)%model )
#                 if defined(compute)
#                   include "compute_coul.f90"
#                 else
#                   include "virial_compute_coul.f90"
#                 endif
                end select
#               if defined(compute)
                  Ecoul = Ecoul + Eij
#               endif
                Wcoul = Wcoul + Wij
                Wsum = Wsum + Wij
              end if
            end associate
            Fij = Wsum*invR2*Rij
            Fi = Fi + Fij
            F(:,j) = F(:,j) - Fij
          end if
        end do
      end associate
      F(:,i) = F(:,i) + Fi
    end do
  end associate
  F = me%Lbox*F
end block
