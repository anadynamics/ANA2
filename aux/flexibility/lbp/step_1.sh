#!/bin/bash

echo -e "Alpha without ligand"
ANA2 avg_4uet.pdb -c config_alfa.cfg -M modes_4uet -Z 5

echo -e "\nAlpha with ligand"
ANA2 avg_4xcp.pdb -c config_alfa.cfg -M modes_4xcp -Z 5

echo -e "\nBeta without ligand"
ANA2 avg_1ifb.pdb -c config_beta.cfg -M modes_1ifb -Z 5

echo -e "\nBeta with ligand"
ANA2 avg_2ifb.pdb -c config_beta.cfg -M modes_2ifb -Z 5

exit 0
