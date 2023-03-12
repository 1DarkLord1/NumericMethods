// Task 3 (Суммирование Кэхена)
#include <vector>
#include <iostream>
#include "utils.hpp"

int main() {
  // Будем сравнивать относительные погрешности naive_sum и kahan_sum на
  // массиве с одинаковыми элементами и фиксированной заранее известной суммой.

  // Возьмём большое n, чтобы в процессе суммирования было задействовано много окон.

  // По результатам видим, что относительна погрешность метода суммирования Кэхана равна 0,
  // тогда как наивное суммирование ошибается на примерно -2% от true_result.

  const float elem = 1.0 / 3;
  const float true_result = 1000000.0;
  const int n = 3000000;
  std::vector<float> x(n, elem);

  float naive_sum_result = utils::naive_sum((const float *) x.data(), n);
  float kahan_sum_result = utils::kahan_sum((const float *) x.data(), n);

  std::cout << "naive_sum relative error: " << utils::relative_error(naive_sum_result, true_result) << std::endl;
  std::cout << "kahan_sum relative error: " << utils::relative_error(kahan_sum_result, true_result) << std::endl;

  return 0;
}