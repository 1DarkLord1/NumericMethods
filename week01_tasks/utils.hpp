// Общие функции, которые используются в разных задачах
#pragma once

namespace utils {
    float naive_sum(const float* x, int n);

    float kahan_sum(const float* x, int n);

    float relative_error(float result, float true_result);
}
