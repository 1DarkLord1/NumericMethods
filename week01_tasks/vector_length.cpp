#include <algorithm>
#include <cmath>
#include <iostream>

namespace vector_length {
    float length(const float *x, int n) {
      if (n == 0) {
        return 0;
      }

      float _max;
      float prev_max;
      float sum = 0;

      for (int i = 0; i < n; i++) {
        if (i == 0) {
          _max = x[i];
          sum = 1;
          continue;
        }

        if (_max < x[i]) {
          prev_max = _max;
          _max = x[i];
          sum *= prev_max / _max;
        }

        sum += (x[i] / _max) * (x[i] / _max);
      }

      return _max * std::sqrt(sum);
    }
}

int main() {
  const int n = 3000000;
  const float elem = 1.0 / 3;

  std::vector<float> x(n, elem);
  std::cout << "vector length: " << vector_length::length((const float *) x.data(), n) << std::endl;
}
