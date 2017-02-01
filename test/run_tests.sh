#!/bin/bash
commands+=("./test_pair_lj_cut small.inp")
commands+=("./test_pair_lj_sf small.inp")
commands+=("./test_pair_softcore_cut small.inp")
commands+=("./test_verlet 2 small.inp")
commands+=("./testfortran 2 small.inp")
commands+=("./test_kspace_ewald water.inp")
#commands+=("")
#commands+=("")
#commands+=("")

for i in "${!commands[@]}"; do
  echo "====================================================================================="
  echo ${commands[$i]}
  echo "-------------------------------------------------------------------------------------"
  eval ${commands[$i]}
  if [ "$?" == "0" ]; then
    echo -e "\n> Passed"
  else
    echo -e "\n> FAILED!"
  fi 
  echo "====================================================================================="
  echo
done

