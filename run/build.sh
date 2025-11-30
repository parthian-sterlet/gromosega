#!/bin/sh

chmod a+x simple.pl
chmod a+x mask.pl
chmod a+x no_mask.pl
chmod a+x command_line_simple
chmod a+x command_line_no_mask
chmod a+x command_line_blacklisted
chmod a+x command_line_whitelisted

g++ -o minimax.exe ../src/minimax.cpp
g++ -o table_rnaseq_filter.exe ../src/table_rnaseq_filter.cpp
g++ -o select_lines01.exe ../src/select_lines01.cpp
g++ -o bed_chr_separation.exe ../src/bed_chr_separation.cpp
g++ -o bed_sort.exe ../src/bed_sort.cpp
g++ -o fasta_muliplefiles.exe ../src/fasta_muliplefiles.cpp
g++ -o fasta_to_plain0.exe ../src/fasta_to_plain0.cpp
g++ -o longext_many.exe ../src/longext_many.cpp
g++ -o mix0.exe ../src/mix0.cpp

chmod a+x minimax.exe
chmod a+x table_rnaseq_filter.exe
chmod a+x select_lines01.exe
chmod a+x bed_chr_separation.exe
chmod a+x bed_sort.exe
chmod a+x fasta_muliplefiles.exe
chmod a+x fasta_to_plain0.exe
chmod a+x longext_many.exe
chmod a+x mix0.exe

