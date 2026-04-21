#!/bin/sh

chmod a+x rnaseq_total_com_line
chmod a+x rnaseq_slice com_line
chmod a+x gmsga_rnaseq_total
chmod a+x gmsga_rnaseq_slice

cd ../src

g++ -o minimax.exe minimax.cpp
g++ -o minimax_pars.exe minimax_pars.cpp
g++ -o minimax_sort.exe minimax_sort.cpp
g++ -o table_rnaseq_adapted.exe table_rnaseq_adapted.cpp
g++ -o table_rnaseq_knife.exe table_rnaseq_knife.cpp
g++ -o select_lines01.exe select_lines01.cpp

chmod a+x minimax.exe
chmod a+x minimax_pars.exe
chmod a+x minimax_sort.exe
chmod a+x table_rnaseq_adapted.exe
chmod a+x table_rnaseq_knife.exe
chmod a+x select_lines01.exe

cd ..
