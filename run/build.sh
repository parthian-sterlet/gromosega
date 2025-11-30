#!/bin/sh

chmod a+x chipseq_com_line
chmod a+x rnaseq_com_line

cd ..
cd src

g++ -o minimax.exe minimax.cpp
g++ -o table_rnaseq_filter.exe table_rnaseq_filter.cpp
g++ -o select_lines01.exe select_lines01.cpp

chmod a+x minimax.exe
chmod a+x table_rnaseq_filter.exe
chmod a+x select_lines01.exe

cd ..
