#include "main.h"

static void init_array(int ni, int nj, int nk, int nl,
                        double *alpha,
                        double *beta,
                        double A[ni][nk],
                        double B[nk][nj],
                        double C[nj][nl],
                        double D[ni][nl]) {

  int i, j;
  *alpha = 1.5;
  *beta = 1.2;
  for (i = 0; i < ni; i++)
    for (j = 0; j < nk; j++)
      A[i][j] = (double) ((i*j+1) % ni) / ni;

  for (i = 0; i < nk; i++)
    for (j = 0; j < nj; j++)
      B[i][j] = (double) (i*(j+1) % nj) / nj;

  for (i = 0; i < nj; i++)
    for (j = 0; j < nl; j++)
      C[i][j] = (double) ((i*(j+3)+1) % nl) / nl;

  for (i = 0; i < ni; i++)
    for (j = 0; j < nl; j++)
      D[i][j] = (double) (i*(j+2) % nk) / nk;
}

static void kernel_2mm( int ni, int nj, int nk, int nl,
                        double alpha,
                        double beta,
                        double tmp[ni][nj],
                        double A[ni][nk],
                        double B[nk][nj],
                        double C[nj][nl],
                        double D[ni][nl]) {

 int i, j, k;
 for (i = 0; i < ni; i++)
        for (j = 0; j < nj; j++) {
                tmp[i][j] = 0.0;
                for (k = 0; k < nk; k++)
                        tmp[i][j] += alpha * A[i][k] * B[k][j];
        }

  for (i = 0; i < ni; i++)
        for (j = 0; j < nl; j++) {
                D[i][j] *= beta;
                for (k = 0; k < nj; k++)
                        D[i][j] += tmp[i][k] * C[k][j];
        }
}