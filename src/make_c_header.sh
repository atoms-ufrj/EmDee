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

function get_parameters {
  grep -i -A100 -e "^\s*type\s*\,\s*extends.*$1" src/$1.f90 |
  grep -i -m1 -B100 "^\s*end\s*type" |
  sed -e "s/\s*//g" |
  grep -i -m1 -e "^real(rb)" |
  sed -e "s/^real(rb)//I" -e "s/:://" -e "s/!.*//" |
  sed -e "s/,/, /g"
}

cat src/emdee_header.h
for model in "$@"; do
    params=$(get_parameters $model | sed -e "s/\([a-zA-Z][a-zA-Z0-9_]*\)/double \1/g" )
    echo "void* EmDee_$model( $params );"
done

