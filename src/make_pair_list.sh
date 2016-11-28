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

# Look for independent modules:
list=""
others=""
search="-v"
for model in "$@"; do
    if [ $(grep -c -e "use\s*pair_" src/${model}.f90) -eq 0 ]; then
        list="$list $model"
        search="$search -e \"use\s*${model}_module\""
    else
        others="$others $model"
    fi
done

# Look for modules that depend solely on those already listed:
until [[ -z $others ]]; do
    group=$others
    others=""
    for model in $group; do
        if [ $(eval "grep $search src/$model.f90" | grep -c -e "use\s*pair_") -eq 0 ]; then
            list="$list $model"
            search="$search -e \"use\s*${model}_module\""
        else
            others="$others $model"
        fi
    done
done

# Output ordered list:
echo $list

