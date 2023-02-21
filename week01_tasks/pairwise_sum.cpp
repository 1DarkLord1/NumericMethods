// Task 04 (Попарное суммирование на SIMD)
#include <cstddef>
#include <vector>
#include <chrono>
#include <iostream>
#include "utils.hpp"

namespace pairwise {

    float pairwise_sum_simd(float *x, int n) {
      for (int len = n; len > 1; len >>= 1) {
        int new_len = len >> 1;

#pragma omp simd
        for (size_t i = 0; i < new_len; i++) {
          x[i] = x[i << 1] + x[(i << 1) + 1];
        }

        if (len & 1) {
          x[new_len - 1] += x[len - 1];
        }
      }
      return x[0];
    }
}

int main() {
  // Будем сравнивать относительные погрешности naive_sum и pairwise_sum_simd на
  // массиве с одинаковыми элементами и фиксированной заранее известной суммой.

  // Возьмём большое n = 30 000 000, чтобы в процессе суммирования было задействовано много окон.

  // По результатам видим, что относительная погрешность метода попарного суммирования равна 0,
  // тогда как наивное суммирование ошибается на примерно -16% от true_result. Оба метода в среднем работают примерно
  // 117 ms (даже если использовать simd в pairwise_sum_simd; видимо это происходит, потому что компилятор
  // оптимизирует цикл в naive_sum).

  const float elem = 1.0 / 3;
  const float true_result = 10000000.0;
  const size_t n = 30000000;
  const int tries = 100;

  std::vector<float> x(n);
  for (size_t i = 0; i < n; i++) {
    x[i] = elem;
  }

  float naive_sum_result = 0;
  float naive_mean_time = 0;

  for (int i = 0; i < tries; i++) {
    std::chrono::steady_clock::time_point naive_begin = std::chrono::steady_clock::now();
    naive_sum_result = utils::naive_sum((const float *) x.data(), n);
    std::chrono::steady_clock::time_point naive_end = std::chrono::steady_clock::now();
    naive_mean_time += (float) std::chrono::duration_cast<std::chrono::milliseconds>(naive_end - naive_begin).count();
  }

  naive_mean_time /= (float) tries;

  float pairwise_sum_result = 0;
  float pairwise_mean_time = 0;

  for (int i = 0; i < tries; i++) {
    std::vector<float> x_ = x;
    std::chrono::steady_clock::time_point pairwise_begin = std::chrono::steady_clock::now();
    pairwise_sum_result = pairwise::pairwise_sum_simd(x_.data(), n);
    std::chrono::steady_clock::time_point pairwise_end = std::chrono::steady_clock::now();
    pairwise_mean_time += (float) std::chrono::duration_cast<std::chrono::milliseconds>(
            pairwise_end - pairwise_begin).count();
  }

  pairwise_mean_time /= (float) tries;

  std::cout << "naive_sum relative error: " << utils::relative_error(naive_sum_result, true_result) <<
            std::endl << "naive_sum mean time elapsed: " << naive_mean_time << "[ms]"
            << std::endl;

  std::cout << "pairwise_sum_simd relative error: " << utils::relative_error(pairwise_sum_result, true_result) <<
            std::endl << "pairwise_sum_simd mean time elapsed: " << pairwise_mean_time << "[ms]"
            << std::endl;

  return 0;
}