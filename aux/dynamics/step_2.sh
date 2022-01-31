#!/bin/bash

#ANA2 4xcp.pdb -d lv_4xcp.nc -c config_2.cfg -t ch_test_lv_2
ANA2 4xcp.pdb -d lv_4xcp.nc -c config_2.cfg -o 2_lv_volume

#ANA2 4xcp.pdb -d hv_4xcp.nc -c config_2.cfg -t ch_test_hv_2
ANA2 4xcp.pdb -d hv_4xcp.nc -c config_2.cfg -o 2_hv_volume

exit 0
