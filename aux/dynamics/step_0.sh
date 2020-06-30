#!/bin/bash

ANA2 4xcp_lv.pdb -c config_1.cfg -f lv_cavity
ANA2 4xcp_hv.pdb -c config_1.cfg -f hv_cavity

ANA2 4xcp_lv.pdb -c config_2.cfg -f lig_lv_cavity
ANA2 4xcp_hv.pdb -c config_2.cfg -f lig_hv_cavity
