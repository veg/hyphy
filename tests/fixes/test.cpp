#include <cstdio>

class _Formula {
  long field;
};

int main(void) {
  _Formula **T = new _Formula *[10]();
  for (int i = 0; i < 10; i++) {
    printf("%x\n", T[i]);
  }
  return 0;
}