#!/bin/bash
nthreads=2
#commands+=("./test_pair_lj_cut $nthreads lj_sample.inp")
#commands+=("./test_pair_lj_sf $nthreads lj_sample.inp")
#commands+=("./test_pair_lj_smoothed $nthreads lj_sample.inp")
#commands+=("./test_pair_lj_shifted_smoothed $nthreads lj_sample.inp")
#commands+=("./test_pair_softcore_cut $nthreads lj_sample.inp")
#commands+=("./test_verlet $nthreads lj_sample.inp")
#commands+=("./test_rigid_body_setup $nthreads spce_sample.inp")
#commands+=("./test_coul_sf $nthreads spce_sample.inp")
#commands+=("./test_coul_damped_smoothed $nthreads spce_sample.inp")
#commands+=("./test_coul_shifted_smoothed $nthreads spce_sample.inp")
#commands+=("./test_coul_long $nthreads spce_sample.inp")
commands+=("./test_rigid_body_exact $nthreads spce_sample.inp")
commands+=("./test_rigid_body_miller $nthreads spce_sample.inp")
#commands+=("./test_respa $nthreads spce_sample.inp")
#commands+=("./test_phase_space_sharing $nthreads spce_sample.inp")
#commands+=("./testfortran $nthreads lj_sample.inp")

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
