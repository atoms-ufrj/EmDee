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

grep -v -e "end\s*interface" src/emdee_header.f03
for model in "$@"; do
    params=$(grep --ignore-case -A100 -e "^\s*type\s*\,\s*extends.*$model" src/$model.f90 | \
             grep --ignore-case -m1 -e "^\s*real\s*(\s*rb\s*)" | \
             sed -e "s/^\s*real\s*(\s*rb\s*)\s*//I" -e "s/::\s*//" -e "s/\s*!.*//")
    echo "  type(c_ptr) function EmDee_$model( $params ) bind(C,name=\"EmDee_$model\")"
    echo "    import :: c_ptr, c_double"
    echo "    real(c_double), value :: $params"
    echo "  end function EmDee_$model"
    echo
done
echo "end interface"
