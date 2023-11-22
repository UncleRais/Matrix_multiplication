#include "main.h"

double bench_t_start, bench_t_end;

static double rtclock() {
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, NULL);
    if (stat != 0)
      printf ("Error return from gettimeofday: %d", stat);
    return (Tp.tv_sec + Tp.tv_usec * 1.0e-6);
}

void bench_timer_start() {
  bench_t_start = omp_get_wtime();
}

void bench_timer_stop() {
  bench_t_end = omp_get_wtime();
}

double bench_timer_get() {
  return(bench_t_end - bench_t_start);
}

double bench_timer_print() {
  printf ("Time in seconds = %0.6lf\n", bench_t_end - bench_t_start);
  return(bench_t_end - bench_t_start);
}