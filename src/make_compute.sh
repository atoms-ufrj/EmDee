#!/bin/bash

#   This file is part of EmDee.
#
#    EmDee is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    EmDee is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with EmDee. If not, see <http://www.gnu.org/licenses/>.
#
#    Author: Charlles R. A. Abreu (abreu@eq.ufrj.br)
#            Applied Thermodynamics and Molecular Simulation
#            Federal University of Rio de Janeiro, Brazil

#---------------------------------------------------------------------------------------------------

if [[ "$1" == "pair" ]] || [[ "$1" == "coul" ]]; then
  echo "type is ($1_none)"
  echo "  Eij = zero"
  echo "  Wij = zero"
  echo
fi

for model in "${@:2}"; do
  echo "type is ($model)"
  echo "  block"
  grep -i -A100 -e "^\s*subroutine\s*${model}_compute" src/$model.f90 |
  grep -i -m1 -B100 "^\s*end\s*subroutine" |
  grep -v -i -e "subroutine" -e "intent(.*)"
  echo "  end block"
  echo
done
