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

echo "module models"
for var in "$@"; do
    echo "  use ${var}_module"
done
for var in pair bond angle dihedral; do
    echo "  use ${var}ModelClass"
done
echo "contains"
echo "!----------"
for var in "$@"; do
    params=$(grep --ignore-case -A100 -e "^\s*type\s*\,\s*extends.*$var" src/$var.f90 | \
    grep --ignore-case -m1 -e "^\s*real\s*(\s*rb\s*)" | \
    sed -e "s/^\s*real\s*(\s*rb\s*)\s*//I" -e "s/::\s*//" -e "s/\s*!.*//")
    echo "  type(c_ptr) function EmDee_$var($params) bind(C,name=\"EmDee_$var\")"
    echo "    real(c_double), value :: $params"
    echo "    type($var), pointer :: model"
    echo "    allocate(model)"
    echo "    call model % setup( [$params] )"
    echo "    EmDee_$var = model % deliver()"
    echo "  end function EmDee_$var"
    echo "!----------"
done
echo "end module models"
