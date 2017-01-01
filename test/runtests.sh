#!/bin/bash
./test_pair_coul_sf ../examples/small.inp
./test_pair_lj_cut_coul_cut ../examples/small.inp
./test_pair_lj_cut_coul_sf ../examples/small.inp
./test_pair_lj_cut ../examples/small.inp
./test_pair_lj_sf_coul_sf ../examples/small.inp
./test_pair_lj_sf ../examples/small.inp
./test_pair_softcore_cut_coul_sf ../examples/small.inp
./test_pair_softcore_cut ../examples/small.inp
./test_pair_softcore_sf_coul_sf ../examples/small.inp
./testfortran 2 ../examples/small.inp

