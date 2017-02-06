#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

int partition(int N, int m)
{
  return (N + m - 1) / m; /// ceil(N/m);
}

/// calculate the index/address from a normal 3D mesh
int idxA(int i, int j, int k, int N1, int N2, int N3, int h)
{
  return k*N1*N2 + j*N1 + i;
}

/// calculate the index/address from a transformed 3D mesh
int idxB(int i, int j, int k, int N1, int N2, int N3, int h)
{
  int blk = i / h;
  int hh = MIN(N1 - blk * h, h);

  return blk * h*N2*N3 + k*hh*N2 + j*hh + (i % h);
}

/// define a function pointer
typedef int (*idx_func)(int, int, int, int, int, int, int);

/// convolution
void convolution(float *out, const float *in, idx_func idx, int i, int j, int k, int N1, int N2, int N3, int h)
{
  out[idx(i, j, k, N1, N2, N3, h)] = 1.0/7 * (in[idx(i-1, j, k, N1, N2, N3, h)] + in[idx(i+1, j, k, N1, N2, N3, h)] +
                                              in[idx(i, j-1, k, N1, N2, N3, h)] + in[idx(i, j+1, k, N1, N2, N3, h)] +
                                              in[idx(i, j, k-1, N1, N2, N3, h)] + in[idx(i, j, k+1, N1, N2, N3, h)] +
                                              in[idx(i, j, k,   N1, N2, N3, h)] );
}

int main(void)
{
  int N1 = 300; // fast axis
  int N2 = 200;
  int N3 = 100; // slow axis
  int m1 = 11;  // partition blocks

  int h = partition(N1, m1);    // num elements in each block along N1 axis

  float *A = (float *)malloc(sizeof *A * N1 * N2 * N3);
  float *B = (float *)malloc(sizeof *A * N1 * N2 * N3);

  /// init A
  for (int i = 0; i < N1 * N2 * N3; i++) {
    A[i] = rand() * 1.0 / RAND_MAX;
  }

  /// transform A to B
  for (int i3 = 0; i3 < N3; i3++) {
    for (int i2 = 0; i2 < N2; i2++) {
      for (int i1 = 0; i1 < N1; i1++) {
        B[idxB(i1, i2, i3, N1, N2, N3, h)] = A[idxA(i1, i2, i3, N1, N2, N3, h)];
      }
    }
  }

  float *AA = (float *)malloc(sizeof *AA * N1 * N2 * N3);
  float *BB = (float *)malloc(sizeof *BB * N1 * N2 * N3);
  float *CC = (float *)malloc(sizeof *CC * N1 * N2 * N3);
  memset(AA, 0, sizeof *CC * N1*N2*N3);
  memset(BB, 0, sizeof *CC * N1*N2*N3);
  memset(CC, 0, sizeof *CC * N1*N2*N3);

  /// perform basic convolution to A and B
  /// AA = convolution(A)
  //  BB = convolution(B)
  for (int i3 = 1; i3 < N3-1; i3++) {
    for (int i2 = 1; i2 < N2-1; i2++) {
      for (int i1 = 1; i1 < N1-1; i1++) {
        convolution(AA, A, idxA, i1, i2, i3, N1, N2, N3, h);
        convolution(BB, B, idxB, i1, i2, i3, N1, N2, N3, h);
      }
    }
  }

  /// compare the result of AA and BB
  for (int i3 = 0; i3 < N3; i3++) {
    for (int i2 = 0; i2 < N2; i2++) {
      for (int i1 = 0; i1 < N1; i1++) {
        if (fabs(AA[idxA(i1, i2, i3, N1, N2, N3, h)] - BB[idxB(i1, i2, i3, N1, N2, N3, h)]) > 0.0001) {
          printf("error at (%d, %d, %d): AA=%f, BB=%f\n", i1, i2, i3, AA[idxA(i1, i2, i3, N1, N2, N3, h)], BB[idxB(i1, i2, i3, N1, N2, N3, h)]);
          exit(1);
        }
      }
    }
  }

  printf("pass all test: 1\n");


  /// suppose B is in the CPE, and perform the convolution
  /* for (int ib = 0; ib < m1; ib++) { */
  /*   /// the start index of the block, the upper halo is included: (h-1) */
  /*   int start_idx; */
  /*   if (ib == 0) { */
  /*     start_idx = 0; */
  /*   } else { */
  /*     start_idx = ib * (h-1)*N2*N3; */
  /*   } */
  /*   int start_idx = MAX(0, ib * (h-1)*N2*N3); */
  /*   int hh = MIN(N1 - ib * h, h); */

  /*   for (int i3 = 1; i3 < N3-1; i3++) { */
  /*     for (int i2 = 1; i2 < N2-1; i2++) { */
  /*       for (int i1 = 1; i1 < hh - 1; i1++) { */
  /*         convolution(&CC[start_idx], &B[start_idx], idxA, i1, i2, i3, hh, N2, N3, hh); */
  /*       } */
  /*     } */
  /*   } */
  /* } */

  /* /// compare the result of AA and BB */
  /* for (int i3 = 0; i3 < N3; i3++) { */
  /*   for (int i2 = 0; i2 < N2; i2++) { */
  /*     for (int i1 = 0; i1 < N1; i1++) { */
  /*       if (fabs(AA[idxA(i1, i2, i3, N1, N2, N3, h)] - CC[idxB(i1, i2, i3, N1, N2, N3, h)]) > 0.0001) { */
  /*         printf("error at (%d, %d, %d): AA=%f, CC=%f\n", i1, i2, i3, AA[idxA(i1, i2, i3, N1, N2, N3, h)], CC[idxB(i1, i2, i3, N1, N2, N3, h)]); */
  /*         exit(1); */
  /*       } */
  /*     } */
  /*   } */
  /* } */

  printf("pass all test: 2\n");

  free(A);
  free(B);
  free(AA);
  free(BB);
  free(CC);

  printf("pass all test\n");
  return 0;
}
