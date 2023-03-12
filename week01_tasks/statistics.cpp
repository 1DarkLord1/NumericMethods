#include <iostream>
#include "statistics.hpp"
#include "utils.hpp"

namespace stat {
    void Statistics::update(float x) {
      _storage.push_back(x);

      _min = _storage.size() == 1 ? x : std::min(_min, x);
      _max = _storage.size() == 1 ? x: std::max(_max, x);

      volatile float y = x - _sum_corr;
      volatile float t = _sum + y;
      _sum_corr = (t - _sum) - y;
      _sum = t;

      volatile float prelast_mean = _mean;
      _mean += (x - _mean) / (float) _storage.size();
      _variance += (x - prelast_mean) * (x - _mean);
    }

    int Statistics::count() const noexcept {
      return (int) _storage.size();
    }

    float Statistics::min() const noexcept {
      return _min;
    }

    float Statistics::max() const noexcept {
      return _max;
    }

    float Statistics::sum() const noexcept {
      return _sum;
    }

    float Statistics::mean() const noexcept {
      return _mean;
    }

    float Statistics::variance() const noexcept {
      return _variance / (float) _storage.size();
    }
}

int main() {
  const int n = 3000000;
  const float elem = 1.0 / 3;

  stat::Statistics s{};
  for (int i = 0; i < n; i++) {
    s.update(elem);
    s.update(2 * elem);
  }

  std::cout << "count: " << s.count() << std::endl <<
            "min: " << s.min() << std::endl <<
            "max: " << s.max() << std::endl <<
            "sum: " << s.sum() << std::endl <<
            "mean: " << s.mean() << std::endl <<
            "variance: " << s.variance() << std::endl;
  return 0;
}