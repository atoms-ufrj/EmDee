#!/bin/bash
commands+=("./test_pair_lj_cut 2 lj_sample.inp")
commands+=("./test_pair_lj_sf 2 lj_sample.inp")
commands+=("./test_pair_softcore_cut 2 lj_sample.inp")
commands+=("./test_verlet 2 lj_sample.inp")
commands+=("./test_rigid_body_setup 2 spce_sample.inp")
commands+=("./test_coul_sf 2 spce_sample.inp")
commands+=("./test_pair_mie_cut 2 lj_sample.inp")
commands+=("./test_coul_long 2 spce_sample.inp")
commands+=("./test_kspace_spme 2 spce_sample.inp")
commands+=("./test_rigid_body_exact 2 spce_sample.inp")
commands+=("./test_rigid_body_miller 2 spce_sample.inp")
#commands+=("./testfortran 2 lj_sample.inp")
#commands+=("./test_kspace_ewald water.inp")
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

