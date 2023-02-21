// Task 05 (Статистика)
#pragma once

#include <vector>

namespace stat {
    class Statistics {
    public:
        void update(float x);            // добавить новый элемент
        int count() const noexcept;

        float min() const noexcept;

        float max() const noexcept;

        float sum() const noexcept;

        float mean() const noexcept;     // среднее
        float variance() const noexcept; // дисперсия
    private:
        std::vector<float> _storage;
        mutable float _min, _max, _sum, _sum_corr, _mean, _variance;
    };
}