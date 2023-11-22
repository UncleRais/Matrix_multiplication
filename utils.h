#include "main.h"

static void print_array(int ni, int nl, double D[ni][nl]) {
  fprintf(stderr, "==BEGIN DUMP_ARRAYS==\n");
  fprintf(stderr, "begin dump: %s", "D");
  fprintf (stderr, "\n");
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
  		fprintf (stderr, "%0.4lf ", D[i][j]);
  	}
  	fprintf (stderr, "\n");
  }
  fprintf(stderr, "\nend   dump: %s\n", "D");
  fprintf(stderr, "==END   DUMP_ARRAYS==\n");
}

static double matrix_equal(int ni, int nj, double tmp1[ni][nj], double tmp2[ni][nj]) {
  int i, j;
  double norm = 0.0;
  for (i = 0; i < ni; ++i) 
    for (j = 0; j < nj; ++j) 
      norm += fabs(tmp1[i][j] - tmp2[i][j]);
  return(norm);
}

static void copy_data(int n, int m, double from[n][m], double to[n][m]) {
  int i, j;
  for (i = 0; i < n; ++i)
    for (j = 0; j < m; ++j) 
      to[i][j] = from[i][j];
}

static void file_printf(FILE* file, char* message, char* str_data, double num_data) {
  if (fabs(num_data + 1.0) < 1e-5) {
  	sprintf(message, "%s %s", str_data, "\n");
  }
  else {
  	sprintf(message, "%s %f %s", str_data, num_data, "\n");
  }
  fputs(message, file);
}

