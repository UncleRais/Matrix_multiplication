#include "initial.h"
#include "initial_optimized.h"
#include "block_multiplication.h"
#include "parallel_for.h"
#include "parallel_task.h"
#include "utils.h"
#include "clock.h"


int main(int argc, char** argv)
{
  int iterations = atoi(argv[1]); 
  int NUM_THREADS = atoi(argv[2]);
  char* filename = argv[3];
  char (*message) [100]; message = (char(*)[100]) malloc(100 * sizeof(char));
  FILE *file = fopen(filename, "at");
  if (!file) {
    printf("Error with file opening. \n");
    return(0);
  }

  int ni = NI;
  int nj = NJ;
  int nk = NK;
  int nl = NL;

  time_t mytime = time(NULL);
  struct tm *now = localtime(&mytime);
  printf("Date: %d.%d.%d\n", now->tm_mday, now->tm_mon + 1, now->tm_year + 1900);
  printf("Time: %d:%d:%d\n", now->tm_hour, now->tm_min, now->tm_sec);

  file_printf(file, *message, "---------------------------------------------------------------", -1);
  sprintf(*message, "Date: %d.%d.%d\n", now->tm_mday, now->tm_mon + 1, now->tm_year + 1900);
  fputs(*message, file);
  sprintf(*message, "Time: %d:%d:%d\n", now->tm_hour, now->tm_min, now->tm_sec, "\n");
  fputs(*message, file);

  file_printf(file, *message, "ni = ", ni);
  file_printf(file, *message, "nj = ", nj);
  file_printf(file, *message, "nl = ", nl);
  file_printf(file, *message, "nk = ", nk);

  double alpha;
  double beta;
  double (*tmp)[ni][nj]; tmp = (double(*)[ni][nj]) malloc((ni) * (nj) * sizeof(double));
  double (*A)[ni][nk];   A = (double(*)[ni][nk]) malloc((ni) * (nk) * sizeof(double));
  double (*B)[nk][nj];   B = (double(*)[nk][nj]) malloc((nk) * (nj) * sizeof(double));
  double (*C)[nj][nl];   C = (double(*)[nj][nl]) malloc((nj) * (nl) * sizeof(double));
  double (*D)[ni][nl];   D = (double(*)[ni][nl]) malloc((ni) * (nl) * sizeof(double));

  double (*B_T)[nj][nk]; B_T = (double(*)[nj][nk]) malloc((nk) * (nj) * sizeof(double));
  double (*C_T)[nl][nj]; C_T = (double(*)[nl][nj]) malloc((nj) * (nl) * sizeof(double));

  double times[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; 

  double (*tmp_reference)[ni][nj]; tmp_reference = (double(*)[ni][nj]) malloc((ni) * (nj) * sizeof(double));
  double (*D_reference)[ni][nl]; D_reference = (double(*)[ni][nl]) malloc((ni) * (nl) * sizeof(double));

  file_printf(file, *message, "Number of iterations = ", iterations);


  int MAX_THREADS = omp_get_max_threads();
  if (MAX_THREADS < NUM_THREADS) 
    omp_set_num_threads(MAX_THREADS);
  else
    omp_set_num_threads(NUM_THREADS);
  printf("MAX_THREADS = %d", MAX_THREADS);
  printf("\n");
  printf("NUM_THREADS = %d", NUM_THREADS);
  printf("\n\n");

  file_printf(file, *message, "MAX_THREADS = ", MAX_THREADS);
  file_printf(file, *message, "NUM_THREADS = ", NUM_THREADS);

//////// SEQUENTIAL
  file_printf(file, *message, "Starting sequential calculation.", -1.0);
  printf("Starting sequential calculation. \n");
  for (int iter = 0; iter < iterations; ++iter) {
    bench_timer_start();
  
    init_array(ni, nj, nk, nl, &alpha, &beta, *A, *B, *C, *D);

    kernel_2mm(ni, nj, nk, nl, alpha, beta, *tmp_reference, *A, *B, *C, *D);
  
    bench_timer_stop();
    times[0] += bench_timer_print();
    file_printf(file, *message, "-iteration time = ", bench_timer_get());
  }
  times[0] = times[0] / iterations;
  file_printf(file, *message, "Sum(times) / iterations = mean_time =", times[0]);
  printf("Sum(times) / iterations = mean_time =  %f", times[0]);
  printf("\n");
  printf("Sequential calculation is finished. \n\n");
////////

  copy_data(ni, nl, *D, *D_reference);

//////// OPTIMIZED
  file_printf(file, *message, "Starting sequential optimized calculation.", -1.0);
  printf("Starting sequential optimized calculation. \n");
  for (int iter = 0; iter < iterations; ++iter) {
    bench_timer_start();
  
    init_array_optimized(ni, nj, nk, nl, &alpha, &beta, *A, *B_T, *C_T, *D);
  
    kernel_2mm_optimized(ni, nj, nk, nl, alpha, beta, *tmp, *A, *B_T, *C_T, *D);
    
    bench_timer_stop();
    times[1] += bench_timer_print();
    file_printf(file, *message, "-iteration time = ", bench_timer_get());

  }
  times[1] = times[1] / iterations;
  file_printf(file, *message, "Sum(times) / iterations = mean_time =", times[1]);
  printf("Sum(times) / iterations = mean_time =  %f", times[1]);
  printf("\n");
  printf("Sequential optimized calculation is finished. \n\n");
  printf("Sum_{i,j}|tmp_optimized[i][j] - tmp_sequential[i][j]| =  %f", matrix_equal(ni, nj, *tmp, *tmp_reference));
  printf("\n");
  printf("Sum_{i,j}|D_optimized[i][j] - D_sequential[i][j]| =  %f", matrix_equal(ni, nl, *D, *D_reference));
  printf("\n\n");
////////

//////// BLOCK SEQUENTIAL
  file_printf(file, *message, "Starting sequential block_mult calculation.", -1.0);
  printf("Starting sequential block_mul calculation. \n");
  for (int iter = 0; iter < iterations; ++iter) {
    bench_timer_start();
  
    init_array_block(ni, nj, nk, nl, &alpha, &beta, *tmp, *A, *B_T, *C_T, *D);
  
    kernel_2mm_block(ni, nj, nk, nl, alpha, beta, *tmp, *A, *B_T, *C_T, *D);
    
    bench_timer_stop();
    times[2] += bench_timer_print();
    file_printf(file, *message, "-iteration time = ", bench_timer_get());

  }
  times[2] = times[2] / iterations;
  file_printf(file, *message, "Sum(times) / iterations = mean_time =", times[2]);
  printf("Sum(times) / iterations = mean_time =  %f", times[2]);
  printf("\n");
  printf("Sequential optimized calculation is finished. \n\n");
  printf("Sum_{i,j}|tmp_block[i][j] - tmp_sequential[i][j]| =  %f", matrix_equal(ni, nj, *tmp, *tmp_reference));
  printf("\n");
  printf("Sum_{i,j}|D_block[i][j] - D_sequential[i][j]| =  %f", matrix_equal(ni, nl, *D, *D_reference));
  printf("\n\n");
////////

//////// OPTMIZED FOR
  file_printf(file, *message, "Starting parallel(for) calculation.", -1.0);
  printf("Starting parallel(for) calculation. \n");
  for (int iter = 0; iter < iterations; ++iter) {
    bench_timer_start();
  
    init_array_for(ni, nj, nk, nl, &alpha, &beta, *A, *B_T, *C_T, *D);
  
    kernel_2mm_for(ni, nj, nk, nl, alpha, beta, *tmp, *A, *B_T, *C_T, *D);

    bench_timer_stop();
    times[3] += bench_timer_print();
    file_printf(file, *message, "-iteration time = ", bench_timer_get());
  }
  times[3] = times[3] / iterations;
  file_printf(file, *message, "Sum(times) / iterations = mean_time =", times[3]);
  printf("Sum(times) / iterations = mean_time =  %f", times[3]);
  printf("\n");
  printf("Parallel(for) calculation is finished. \n\n");
  printf("Sum_{i,j}|tmp_for[i][j] - tmp_sequential[i][j]| =  %f", matrix_equal(ni, nj, *tmp, *tmp_reference));
  printf("\n");
  printf("Sum_{i,j}|D_for[i][j] - D_sequential[i][j]| =  %f", matrix_equal(ni, nl, *D, *D_reference));
  printf("\n\n");
////////

//////// OPTIMIZED TASK
  file_printf(file, *message, "Starting parallel(task) calculation.", -1.0);
  printf("Starting parallel(task) calculation. \n");
  for (int iter = 0; iter < iterations; ++iter) {
    bench_timer_start();
    
    init_array_task(ni, nj, nk, nl, &alpha, &beta, *A, *B_T, *C_T, *D);
  
    kernel_2mm_task(ni, nj, nk, nl, alpha, beta, *tmp, *A, *B_T, *C_T, *D);
                   
    bench_timer_stop();
    times[4] += bench_timer_print();
    file_printf(file, *message, "-iteration time = ", bench_timer_get());
  }
  times[4] = times[4] / iterations;
  file_printf(file, *message, "Sum(times) / iterations = mean_time =", times[4]);
  printf("Sum(times) / iterations = mean_time =  %f", times[4]);
  printf("\n");
  printf("Parallel(task) calculation is finished. \n\n");

  printf("Sum_{i,j}|tmp_task[i][j] - tmp_sequential[i][j]| =  %f", matrix_equal(ni, nj, *tmp, *tmp_reference));
  printf("\n");
  printf("Sum_{i,j}|D_task[i][j] - D_sequential[i][j]| =  %f", matrix_equal(ni, nl, *D, *D_reference));
  printf("\n\n");
////////

//////// BLOCK FOR
  file_printf(file, *message, "Starting parallel_block(for) calculation.", -1.0);
  printf("Starting parallel_block(for) calculation. \n");
  for (int iter = 0; iter < iterations; ++iter) {
    bench_timer_start();
    
    init_array_block_for(ni, nj, nk, nl, &alpha, &beta, *tmp, *A, *B_T, *C_T, *D);
  
    kernel_2mm_block_for(ni, nj, nk, nl, alpha, beta, *tmp, *A, *B_T, *C_T, *D);
                   
    bench_timer_stop();
    times[5] += bench_timer_print();
    file_printf(file, *message, "-iteration time = ", bench_timer_get());
  }
  times[5] = times[5] / iterations;
  file_printf(file, *message, "Sum(times) / iterations = mean_time =", times[5]);
  printf("Sum(times) / iterations = mean_time =  %f", times[5]);
  printf("\n");
  printf("Parallel_block(for) calculation is finished. \n\n");

  printf("Sum_{i,j}|tmp_block_for[i][j] - tmp_sequential[i][j]| =  %f", matrix_equal(ni, nj, *tmp, *tmp_reference));
  printf("\n");
  printf("Sum_{i,j}|D_block_for[i][j] - D_sequential[i][j]| =  %f", matrix_equal(ni, nl, *D, *D_reference));
  printf("\n\n");
////////

//////// BLOCK TASK
  file_printf(file, *message, "Starting parallel_block(task) calculation.", -1.0);
  printf("Starting parallel_block(task) calculation. \n");
  for (int iter = 0; iter < iterations; ++iter) {
    bench_timer_start();
    
    init_array_block_task(ni, nj, nk, nl, &alpha, &beta, *tmp, *A, *B_T, *C_T, *D);
  
    kernel_2mm_block_task(ni, nj, nk, nl, alpha, beta, *tmp, *A, *B_T, *C_T, *D);
                   
    bench_timer_stop();
    times[6] += bench_timer_print();
    file_printf(file, *message, "-iteration time = ", bench_timer_get());
  }
  times[6] = times[6] / iterations;
  file_printf(file, *message, "Sum(times) / iterations = mean_time =", times[6]);
  printf("Sum(times) / iterations = mean_time =  %f", times[6]);
  printf("\n");
  printf("Parallel_block(task) calculation is finished. \n\n");

  printf("Sum_{i,j}|tmp_block_task[i][j] - tmp_sequential[i][j]| =  %f", matrix_equal(ni, nj, *tmp, *tmp_reference));
  printf("\n");
  printf("Sum_{i,j}|D_block_task[i][j] - D_sequential[i][j]| =  %f", matrix_equal(ni, nl, *D, *D_reference));
  printf("\n\n");
////////


  file_printf(file, *message, "Initial / Optimized: ", times[0] / times[1]);
  file_printf(file, *message, "Initial / Block: ", times[0] / times[2]);
  file_printf(file, *message, "Optimized / Parallel(for): ", times[1] / times[3]);
  file_printf(file, *message, "Optimized / Parallel(task): ", times[1] / times[4]);
  file_printf(file, *message, "Block / Parallel_block(for): ", times[2] / times[5]);
  file_printf(file, *message, "Block / Parallel_block(task): ", times[2] / times[6]);

  printf("Initial / Optimized: %f", times[0] / times[1]);
  printf("\n");
  printf("Initial / Block: %f", times[0] / times[2]);
  printf("\n");
  printf("Optimized / Parallel(for) : %f", times[1] / times[3]);
  printf("\n");
  printf("Optimized / Parallel(task): %f", times[1] / times[4]);
  printf("\n");
  printf("Block / Parallel_block(for): %f", times[2] / times[5]);
  printf("\n");
  printf("Block / Parallel_block(task): %f", times[2] / times[6]);
  printf("\n\n");

  free((void*)tmp);
  free((void*)A);
  free((void*)B);
  free((void*)C);
  free((void*)D);
  free((void*)tmp_reference);
  free((void*)D_reference);
  free((void*)B_T);
  free((void*)C_T);
  free((void*)message);

  fclose(file);

  return 0;
}
