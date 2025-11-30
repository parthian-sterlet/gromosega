#!/bin/sh

chmod a+x command_line_chipseq
chmod a+x command_line_rnaseq

g++ -o minimax.exe ../src/minimax.cpp
g++ -o table_rnaseq_filter.exe ../src/table_rnaseq_filter.cpp
g++ -o select_lines01.exe ../src/select_lines01.cpp

chmod a+x minimax.exe
chmod a+x table_rnaseq_filter.exe
chmod a+x select_lines01.exe
