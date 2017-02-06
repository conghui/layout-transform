#include <sys/time.h>
#include <math.h>
#include <vector>
#include <iostream>


void compute(float &x) {
    for (int i = 1; i < 10; i++) {
      x += sqrt(sin(x) + x*x / (x+i) + log(x+i) + exp(sqrt(x*(x+i)))) / sqrt(sin(x) + x*x / (x+i) + log(x+i) + exp(sqrt(x*(x+i))));
    }
}

int main(void)
{
  int N = 120;
  int M = 100000;
  std::vector<float> input(N * M);

  /// init
  for (size_t i = 0; i < input.size(); i++) {
    input[i] = 10.0 * rand() / RAND_MAX + 1;
  }

  struct timeval start, stop;
  gettimeofday(&start, NULL);

  /// compute
#pragma omp parallel for
  for (int i = 0; i < N; i++) {
    float *p = &input[i*M];

    /// add your MIC directives here
    for (int j = 0; j < M; j++) {
      compute(p[i]);
    }

  }

  gettimeofday(&stop, NULL);
  float elapse = stop.tv_sec - start.tv_sec + (stop.tv_usec - start.tv_usec)*1e-6;
  std::cout << "time used for two-level loop computation is: " << elapse << std::endl;

  return 0;
}
