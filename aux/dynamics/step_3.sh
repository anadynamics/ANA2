#!/bin/bash

#ANA2 4xcp.pdb -d lv_4xcp.nc -c config_3.cfg -t ch_test_lv_3
ANA2 4xcp.pdb -d lv_4xcp.nc -c config_3.cfg -o 3_lv_volume

#ANA2 4xcp.pdb -d hv_4xcp.nc -c config_3.cfg -t ch_test_hv_3
ANA2 4xcp.pdb -d hv_4xcp.nc -c config_3.cfg -o 3_hv_volume

exit 0
