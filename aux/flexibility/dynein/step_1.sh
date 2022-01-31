#!/bin/bash

echo -e "LC8 cavity C"
ANA2 lc8.pdb -c ecf.cfg -M modes_lc8 -F frequencies_lc8 -Z 5

echo -e "\nLC8 cavity D"
ANA2 lc8.pdb -c edf.cfg -M modes_lc8 -F frequencies_lc8 -Z 5

echo -e "\nTcTex cavity C"
ANA2 tctex.pdb -c acb.cfg -M modes_tctex -F frequencies_tctex -Z 5

echo -e "\nTcTex cavity D"
ANA2 tctex.pdb -c adb.cfg -M modes_tctex -F frequencies_tctex -Z 5

exit 0
