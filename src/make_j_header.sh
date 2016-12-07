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

# Remove end of module:
line=$(grep -n "end" src/julia_wrapper.jl | tail -n1 | cut -f1 -d:)
sed "${line}d" src/julia_wrapper.jl

for model in "$@"; do

  # Read parameter definitions from model file:
  IFS=$'\n' definitions=($(get_definitions $model))

  # Parse parameter definitions:
  params=()
  types=()
  jtypes=()
  array=()
  descriptions=()
  for def in "${definitions[@]}"; do
    params+=($(echo $def | sed -e "s/\s*!.*//" | grep -oE '[^ ]+$' | sed -e "s/(.*)//"))
    type=$(echo $def | grep -i -oE -e "^(real|integer)" | tr '[:upper:]' '[:lower:]')
    types+=($type)
    if [ -z $(echo $def | grep -i -oE "allocatable") ]; then
      array+=(0)
      jtypes+=($(echo $type | sed -e 's/real/Float64/' -e 's/integer/Int/'))
    else
      array+=(1)
      jtypes+=($(echo $type | sed -e 's/real/Vector{Float64}/' -e 's/integer/Vector{Int}/'))
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

  # Write down interface to EmDee_model routine:
  allparams=()
  for ((i=0; i< ${#params[@]}; ++i)); do
    allparams+=("${params[i]}::${jtypes[i]}")
    [[ ${types[i]} == "integer" ]] && [ ${array[i]} -eq 1 ] && params[i]="Vector{Int32}(${params[i]})"
  done
  echo "function $model( $(csv ${allparams[*]}) )"
  if [[ -z ${params[@]} ]]; then
    echo "  return ccall( (:EmDee_$model,\"libemdee\"), Ptr{Void}, () )"
  else
    echo "  return ccall( (:EmDee_$model,\"libemdee\"), Ptr{Void},"
    echo "                ($(csv ${jtypes[*]} | sed -e "s/Vector/Ptr/g" -e "s/Int/Int32/g")),"
    echo "                $(csv ${params[*]}) )"
  fi
  echo "end"
  printf "#%99s\n" | tr ' ' '-'
done

echo "end"
