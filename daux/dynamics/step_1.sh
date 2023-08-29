#!/bin/bash

#ANA2 4xcp.pdb -d lv_4xcp.nc -c config_1.cfg -t ch_test_lv_1
ANA2 4xcp.pdb -d lv_4xcp.nc -c config_1.cfg -o 1_lv_volume

#ANA2 4xcp.pdb -d hv_4xcp.nc -c config_1.cfg -t ch_test_hv_1
ANA2 4xcp.pdb -d hv_4xcp.nc -c config_1.cfg -o 1_hv_volume

exit 0
