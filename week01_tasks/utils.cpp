#include <cstddef>
#include "utils.hpp"

namespace utils {
    float naive_sum(const float* x, int n) {
      volatile float s = 0;

      for (int i = 0; i < n; i++) {
        s += x[i];
      }

      return s;
    }

    float kahan_sum(const float *x, int n) {
      volatile float s = 0;
      volatile float c = 0;

      for (int i = 0; i < n; i++) {
        volatile float y = x[i] - c;
        volatile float t = s + y;
        c = (t - s) - y;
        s = t;
      }

      return s;
    }

    float relative_error(float result, float true_result) {
      return (result - true_result) / true_result;
    }
}