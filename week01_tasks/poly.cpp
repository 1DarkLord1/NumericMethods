// Task 2 (Вычисление полиномов)
#include <iostream>
#include <cmath>

namespace poly {
    float polynomial(float x, const float* a, int n) {
      float result = a[n];

      for (int i = n; i >= 1; i--) {
        result = std::fma(x, result, a[i - 1]);
      }

      return result;
    }
}

int main() {
  const float a[] = {1, 1.5, 3.25, 6};
  std::cout << poly::polynomial(3, a, 3) << std::endl;

  return 0;
}