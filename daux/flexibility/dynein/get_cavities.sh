#!/bin/bash

ANA2 lc8.pdb -c ecf.cfg -f ecf
ANA2 lc8.pdb -c edf.cfg -f edf
ANA2 tctex.pdb -c acb.cfg -f acb
ANA2 tctex.pdb -c adb.cfg -f adb

exit 0
