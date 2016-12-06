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

function get_definitions {
  grep -i -A100 -e "^\s*type\s*\,\s*extends.*$1" src/$1.f90 |
  grep -i -m1 -B100 "^\s*end\s*type" |
  grep -i "!<>" |
  sed "s/^\s*//"
}

echo "module models"
for model in "$@"; do
    echo "use ${model}_module"
done
echo "use dihedralModelClass" # REMOVE WHEN DIHEDRAL MODELS HAVE BEEN ADDED
echo "contains"
echo

for model in "$@"; do

  # Read parameter definitions from model file:
  IFS=$'\n' definitions=($(get_definitions $model))

  # Parse parameter definitions:
  params=()
  types=()
  ftypes=()
  array=()
  descriptions=()
  for def in "${definitions[@]}"; do
    params+=($(echo $def | sed -e "s/\s*!.*//" | grep -oE '[^ ]+$' | sed -e "s/(.*)//"))
    type=$(echo $def | grep -i -oE -e "^(real|integer)" | tr '[:upper:]' '[:lower:]')
    types+=($type)
    ftypes+=($(echo $type | sed -e 's/real/real(c_double)/' -e 's/integer/integer(c_int)/'))
    if [ -z $(echo $def | grep -i -oE "allocatable") ]; then
      array+=(0)
    else
      array+=(1)
    fi
    descriptions+=($(echo $def | grep -oE "!<>.*" | sed "s/!<>\s*//" ))
  done

  # Check if allocatable parameters and their sizes have been properly defined:
  last_type="real"
  last_array="0"
  for ((i=0; i<${#params[@]}; i++)); do
    if [ ${array[i]} -eq 1 ]; then
      if [ $last_array -eq 0 ] && [[ $last_type != "integer" ]]; then
        echo "ERROR: parsing of model $model failed.">&2 && exit 1
      fi
    fi
    last_type=${types[i]}
    last_array=${array[i]}
  done

  # Write down EmDee_model routine:
  allparams=${params[@]}
  echo "  type(c_ptr) function EmDee_$model( ${allparams// /, } ) &"
  echo "    bind(C,name=\"EmDee_$model\")"

  if [[ -z $allparams ]]; then

    echo "    type($model), pointer :: model"
    echo "    allocate(model)"
    echo "    call model % setup()"

  else

    # Declare scalar arguments as real or integer and array arguments as type(c_ptr):
    for ((i=0; i<${#params[@]}; i++)); do
      if [ ${array[i]} -eq 0 ]; then
        echo "    ${ftypes[i]}, value :: ${params[i]}"
      else
        echo "    type(c_ptr), value :: ${params[i]}"
      fi
    done

    # Declare a scalar pointer to model and array pointers to array arguments:
    echo "    type($model), pointer :: model"
    for ((i=0; i<${#params[@]}; i++)); do
      [ ${array[i]} -eq 1 ] && echo "    ${ftypes[i]}, pointer :: ptr_${params[i]}(:)"
    done

    # Allocate model pointer:
    echo "    allocate(model)"
    for ((i=0; i<${#params[@]}; i++)); do
      [ ${array[i]} -eq 1 ] && echo "    call c_f_pointer( ${params[i]}, ptr_${params[i]}, [$size] )"
      [ ${array[i]} -eq 0 ] && [[ ${types[i]} == "integer" ]] && size=${params[i]}
      [ ${array[i]} -eq 1 ] && params[i]="ptr_${params[i]}"
    done

    # Split parameters into real and integer ones:
    intpar=()
    realpar=()
    for ((i=0; i<${#params[@]}; i++)); do
      if [[ ${types[i]} == "integer" ]]; then
        intpar+=(${params[i]})
      else
        realpar+=(${params[i]})
      fi
    done

    # Call model setup routine:
    realparams=${realpar[@]}
    intparams=${intpar[@]}
    if [ -z $intparams ]; then
      echo "    call model % setup( [${realparams// /, }] )"
    elif [ -z $realparams ]; then
      echo "    call model % setup( iparams = [${intparams// /, }] )"
    else
      echo "    call model % setup( [${realparams// /, }], [${intparams// /, }] )"
    fi

  fi

  echo "    EmDee_$model = model % deliver()"
  echo "  end function EmDee_$model"
  echo
done

echo "end module models"
