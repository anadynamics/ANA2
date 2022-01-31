#!/bin/bash

echo -e "LC8 cavity C"
ANA2 lc8.pdb -c ecf.cfg -M modes_lc8 -F frequencies_lc8

echo -e "\nLC8 cavity D"
ANA2 lc8.pdb -c edf.cfg -M modes_lc8 -F frequencies_lc8

echo -e "\nTcTex cavity C"
ANA2 tctex.pdb -c acb.cfg -M modes_tctex -F frequencies_tctex

echo -e "\nTcTex cavity D"
ANA2 tctex.pdb -c adb.cfg -M modes_tctex -F frequencies_tctex

exit 0
