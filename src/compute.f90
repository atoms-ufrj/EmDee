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
  real(rb) :: r2, invR2, invR, Wij, Qi, QiQj, WCij, rFc, Eij, r2fac, u, u2, u3, G, WG
  real(rb) :: Rij(3), Ri(3), Fi(3), Fij(3)
  logical  :: icharged, ijcharged
#if defined(compute)
  integer  :: l, layer
  real(rb) :: ECij
#endif
  real(rb), allocatable :: Rvec(:,:)

  F = zero
#ifndef fast
  Wpair = zero
  Wcoul = zero
#endif
#if defined(compute)
  associate ( neighbor => me%neighbor(thread), &
              LayerEpair => me%threadEpair(:,thread), &
              LayerEcoul => me%threadEcoul(:,thread) )
    Epair = zero
    Ecoul = zero
    LayerEpair = zero
    LayerEcoul = zero
#else
  associate ( neighbor => me%neighbor(thread) )
#endif
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
#if defined(fast)
                  jlist => neighbor%item(neighbor%first(i):neighbor%middle(i)) )
#else
                  jlist => neighbor%item(neighbor%first(i):neighbor%last(i)) )
#endif 
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
#if defined(compute)
                  include "compute_pair.f90"
#else
                  include "virial_compute_pair.f90"
#endif
                end select
                select case (model%modifier)
#if defined(compute)
                  case (1)
                    Eij = Eij + model%eshift
#endif
                  case (2)
                    rFc = model%fshift/invR
                    Wij = Wij - rFc
#if defined(compute)
                    Eij = Eij + model%eshift + rFc
#endif
                  case (3,4)
                    r2fac = r2*L2*model%factor
                    if (r2fac > model%Rm2fac) then
#ifndef compute
                      select type ( model )
                        include "energy_compute_pair.f90"
                      end select
#endif
                      u = r2fac - model%Rm2fac
                      u2 = u*u
                      u3 = u*u2
                      G = 1.0_rb + u3*(15.0_rb*u - 6.0_rb*u2 - 10.0_rb)
                      WG = -60.0_rb*u2*(2.0_rb*u - u2 - 1.0_rb)*r2fac
                      Eij = Eij + model%eshift
                      Wij = Wij*G + Eij*WG
                      Eij = Eij*G
                    end if
                end select
              end associate
#if defined(compute)
              Epair = Epair + Eij
#endif
#ifndef fast
              Wpair = Wpair + Wij
#endif
              if (ijcharged.and.pair%coulomb) then
                QiQj = pair%kCoul*Qi*me%charge(j)
#if defined(compute)
                select type ( model => me%coul(me%layer)%model )
                  include "compute_coul.f90"
                end select
                Ecoul = Ecoul + ECij
#else
                select type ( model => me%coul(me%layer)%model )
                  include "virial_compute_coul.f90"
                end select
#endif
#ifndef fast
                Wcoul = Wcoul + WCij
#endif
                Wij = Wij + WCij
              end if
            end associate
#if defined(fast)
            if (r2 > InRc2) then
              u = me%invDeltaR/invR + me%mInRcByDeltaR
              Wij = Wij*(one + u*u*(two*u - three))
            end if
#endif
            Fij = Wij*invR2*Rij
            Fi = Fi + Fij
            F(:,j) = F(:,j) - Fij
#if defined(compute)
            if (me%multilayer(jtype,itype)) then
              LayerEpair(me%layer) = LayerEpair(me%layer) + Eij
              LayerEcoul(me%layer) = LayerEcoul(me%layer) + ECij
              do l = 1, me%nlayers-1
                layer = me%other_layer(l)
                associate( pair => me%pair(itype,jtype,layer) )
                  select type ( model => pair%model )
                    include "energy_compute_pair.f90"
                  end select
! MODEL MODIFICATION MUST COME HERE
                  if (ijcharged.and.pair%coulomb) then
                    QiQj = pair%kCoul*Qi*me%charge(j)
                    select type ( model => me%coul(me%layer)%model )
                      include "energy_compute_coul.f90"
                    end select
                  end if
                  LayerEpair(layer) = LayerEpair(layer) + Eij
                  LayerEcoul(layer) = LayerEcoul(layer) + ECij
                end associate
              end do
            end if
#endif
          end if
        end do
      end associate
      F(:,i) = F(:,i) + Fi
    end do
  end associate
  F = me%Lbox*F
end block
