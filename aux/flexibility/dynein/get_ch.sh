#!/bin/bash

ANA2 lc8.pdb -c ecf.cfg -t ch_ecf
ANA2 lc8.pdb -c edf.cfg -t ch_edf
ANA2 tctex.pdb -c acb.cfg -t ch_acb
ANA2 tctex.pdb -c adb.cfg -t ch_adb

exit 0
