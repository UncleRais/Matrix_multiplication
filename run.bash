#!/bin/bash
#sed -i 's/\r$//' run.bash
echo "---------- -O3"
gcc -fopenmp -std=c99 -o run_gcc_O3 -O3 "main_kernels.c"
./run_gcc_O3.exe 5 2 "run_gcc_O3_2.txt"
./run_gcc_O3.exe 5 4 "run_gcc_O3_4.txt"
./run_gcc_O3.exe 5 6 "run_gcc_O3_6.txt"
./run_gcc_O3.exe 5 8 "run_gcc_O3_8.txt"
./run_gcc_O3.exe 5 10 "run_gcc_O3_10.txt"