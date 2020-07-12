#!/bin/bash

ANA2 avg_4uet.pdb -c config_alfa.cfg -t ch_4uet

ANA2 avg_4xcp.pdb -c config_alfa.cfg -t ch_4xcp

ANA2 avg_1ifb.pdb -c config_beta.cfg -t ch_1ifb

ANA2 avg_2ifb.pdb -c config_beta.cfg -t ch_2ifb

exit 0
