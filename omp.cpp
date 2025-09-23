#include <omp.h>
#include <stdio.h>
int main() {
#pragma omp parallel for default(none)
  for (int i = 0; i < 20; i++) {
    int tid = omp_get_thread_num();
    printf("Hello world from omp thread %d\n", tid);
  }

  return 0;
}
