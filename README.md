# Matrix_multiplication

Contains sequential and parallel implementations of following algorithms:

* General matrix multiplication
* Block matrix multiplication

Note that present implementations are intended for academic purposes, as such they are not meant to be used in any sort of high-performance production code.

## Compilation

* Recommended compilers: gcc, icc, clang
* Requires -std=c99
* Requires openmp option -fopenmp or -openmp

* ## Usage

Solving task:<br>
$$T = \alpha \cdot A \cdot B, \quad D = \beta \cdot D + T \cdot C, \quad \alpha, \beta \in \mathbb{R}, $$
where<br>
$T_{n_i \times n_j} , A_{n_i \times n_k} , B_{n_k \times n_j} , C_{n_j \times n_l}, D_{n_i \times n_l}$ - real matrices, $n_i , n_j , n_l , n_k \in \mathbb{N}$.<br>

## Example

Compilation of "main_kernels.c" with level of optimization -O3. The result of code's work is the mean time (sum / number_of_iterations) of algorithms' work. Actually, was used to explore the relationship between the scale of task and number of simulation threads.
```bash
gcc -fopenmp -std=c99 -o run_gcc_O3 -O3 "main_kernels.c"
### <way_to_exe> <number_of_iterations> <number_of_threads> <name_of_output_file> ###
./run_gcc_O3.exe 5 2 "run_gcc_O3_2.txt"
./run_gcc_O3.exe 5 4 "run_gcc_O3_4.txt"
./run_gcc_O3.exe 5 6 "run_gcc_O3_6.txt"
./run_gcc_O3.exe 5 8 "run_gcc_O3_8.txt"
./run_gcc_O3.exe 5 10 "run_gcc_O3_10.txt"
```

Several examples of calling matrix multiplication functions can be found in "main.h" and "main_kernels.h".   
