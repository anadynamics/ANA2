#!/bin/bash

ANA2 input_pdb.pdb -f 5_cavity -d input_trajectory.nc -c config_4.cfg -t convex_hull_testing_5
ANA2 input_pdb.pdb -f 5_cavity -d input_trajectory.nc -c config_4.cfg

exit 0
