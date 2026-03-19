#!/bin/bash
gnuplot gnufile_TTU_GWU.gplot
gnuplot gnufile_TTU_GWU_para3-4.gplot

evince parameter_match_1.eps&
sleep 1
evince parameter_match_2.eps&
sleep 1
evince parameter_match_3.eps&
sleep 1
evince parameter_match_4.eps&
