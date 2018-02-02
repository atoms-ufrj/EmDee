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

function get_definitions {
  grep -i -A100 -e "^\s*type\s*\,\s*extends.*$1" src/$1.f90 |
  grep -i -m1 -B100 "^\s*end\s*type" |
  grep -i "!<>" |
  sed "s/^\s*//"
}

function csv {
  local IFS=","
  echo "$*" | sed "s/,/, /g"
}

echo "#ifdef __cplusplus"
echo "#include <stdbool.h>"
echo "extern \"C\" {"
echo "#endif"
echo

cat src/emdee_header.h

for model in "$@"; do

  # Read parameter definitions from model file:
  IFS=$'\n' definitions=($(get_definitions $model))

  # Parse parameter definitions:
  params=()
  types=()
  ctypes=()
  array=()
  descriptions=()
  for def in "${definitions[@]}"; do
    params+=($(echo $def | sed -e "s/\s*!.*//" | grep -oE '[^ ]+$' | sed -e "s/(.*)//"))
    type=$(echo $def | grep -i -oE -e "^(real|integer)" | tr '[:upper:]' '[:lower:]')
    types+=($type)
    ctypes+=($(echo $type | sed -e 's/real/double/' -e 's/integer/int/'))
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

  allparams=()
  for ((i=0; i< ${#params[@]}; ++i)); do
    if [ ${array[i]} -eq 0 ]; then
      allparams+=("${ctypes[i]} ${params[i]}")
    else
      allparams+=("${ctypes[i]}* ${params[i]}")
    fi
  done
  echo "void* EmDee_$model( $(csv ${allparams[*]}) );"
  echo
done


echo
echo "#ifdef __cplusplus"
echo "}"
echo "#endif"
