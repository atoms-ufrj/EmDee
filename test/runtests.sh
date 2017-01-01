#!/bin/bash
commands+=("./test_pair_coul_sf ../examples/small.inp")
commands+=("./test_pair_lj_cut_coul_cut ../examples/small.inp")
commands+=("./test_pair_coul_sf ../examples/small.inp")
commands+=("./test_pair_lj_cut_coul_cut ../examples/small.inp")
commands+=("./test_pair_lj_cut_coul_sf ../examples/small.inp")
commands+=("./test_pair_lj_cut ../examples/small.inp")
commands+=("./test_pair_lj_sf_coul_sf ../examples/small.inp")
commands+=("./test_pair_lj_sf ../examples/small.inp")
commands+=("./test_pair_softcore_cut_coul_sf ../examples/small.inp")
commands+=("./test_pair_softcore_cut ../examples/small.inp")
commands+=("./test_pair_softcore_sf_coul_sf ../examples/small.inp")
commands+=("./test_verlet 2 ../examples/small.inp")
commands+=("./testfortran 2 ../examples/small.inp")
#commands+=("")
#commands+=("")
#commands+=("")


for i in "${!commands[@]}"; do
  echo "======================================================================"
  echo ${commands[$i]}
  echo "----------------------------------------------------------------------"
  eval ${commands[$i]}
  echo "======================================================================"
  echo
done
