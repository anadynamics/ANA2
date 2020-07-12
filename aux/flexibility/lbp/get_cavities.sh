#!/bin/bash

echo -e "Alpha without ligand"
ANA2 avg_4uet.pdb -c config_alfa.cfg -f cav_4uet

echo -e "\nAlpha with ligand"
ANA2 avg_4xcp.pdb -c config_alfa.cfg -f cav_4xcp

echo -e "\nBeta without ligand"
ANA2 avg_1ifb.pdb -c config_beta.cfg -f cav_1ifb

echo -e "\nBeta with ligand"
ANA2 avg_2ifb.pdb -c config_beta.cfg -f cav_2ifb

exit 0
