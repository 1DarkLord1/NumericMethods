// Tasks 1, 2

#include <random>
#include <functional>
#include <cmath>
#include <iostream>
#include <set>
#include <cassert>

namespace global_opt {
    int iters_amount = 0;

    double
    random_search_opt(const std::function<double(double)> &f, double left, double right, double min_threshold,
                      int max_iters) {
        iters_amount = 0;
        std::mt19937 gen(42);
        std::uniform_real_distribution<double> distr(-1, 1);
        double min_pt = (left + right) / 2;

        for (int i = 0; i < max_iters; i++) {
            if (f(min_pt) < min_threshold) {
                break;
            }

            double delta = distr(gen);
            double pt = min_pt + delta;

            if (pt < left || pt > right) {
                continue;
            }

            double pt_value = f(pt);

            if (pt_value < f(min_pt)) {
                min_pt = pt;
            }

            iters_amount++;
        }

        return min_pt;
    }

    std::vector<double> get_peaks(const std::function<double(double)> &f, double L, const std::set<double>& pts) {
        std::vector<double> peaks;
        auto prev = pts.begin();
        auto next = pts.begin();
        next++;

        while (next != pts.end()) {
            double x = *prev;
            double y = *next;
            peaks.push_back(0.5 / L * (f(y) - f(x)) + 0.5 * (y + x));
            prev++;
            next++;
        }

        return peaks;
    }

    double
    lipo_opt(const std::function<double(double)> &f, double left, double right, double L, double eps_a, double eps_f,
             int max_iters) {
        iters_amount = 0;
        std::mt19937 gen(42);

        std::set<double> pts;
        pts.insert(left);
        pts.insert((left + right) / 2);
        double max_value_pt = *std::max_element(pts.begin(), pts.end(), [&f](double x, double y) {
            return f(x) < f(y);
        });
        double max_value = f(max_value_pt);

        for (int i = 0; i < max_iters; i++) {
            std::vector<double> peaks = get_peaks(f, L, pts);
            std::uniform_int_distribution<int> distr(0, peaks.size() - 1);
            double pt = peaks[distr(gen)];
            double pt_value = f(pt);

            assert(pt >= left || pt <= right);

            double min_sum_pt = *std::min_element(pts.begin(), pts.end(), [&f, L, pt](double x, double y) {
                return f(x) + L * fabs(pt - x) < f(y) + L * fabs(pt - y);
            });

            double min_sum = f(min_sum_pt) + L * fabs(pt - min_sum_pt);

            if (fabs(pt - max_value_pt) < eps_a || fabs(pt_value - max_value) < eps_f) {
                break;
            }

            if (min_sum >= max_value) {
                pts.insert(pt);
                if (max_value < pt_value) {
                    max_value = pt_value;
                    max_value_pt = pt;
                }
            }
            else {
                i--;
                continue;
            }

            iters_amount++;
        }

        return max_value_pt;
    }

    uint64_t number_to_gray(uint64_t n) {
        return n ^ (n >> 1);
    }

    uint64_t gray_to_number(uint64_t g) {
        uint64_t n = 0;
        for (; g != 0; g >>= 1) { n ^= g; }
        return n;
    }

    uint64_t double_to_bitcode(double x) {
        double *x_ptr = &x;
        return number_to_gray(*reinterpret_cast<uint64_t *>(x_ptr));
    }

    double bitcode_to_double(uint64_t code) {
        uint64_t x = gray_to_number(code);
        uint64_t *x_ptr = &x;
        return *reinterpret_cast<double *>(x_ptr);
    }

    std::pair<uint64_t, uint64_t>
    cross_and_mutate(uint64_t x, uint64_t y, std::uniform_int_distribution<int> &bit_distr,
                     std::mt19937 &gen) {
        int boundary_bit = bit_distr(gen);
        uint64_t yx = y & ((uint64_t) (-1) - ((1 << boundary_bit) - 1)) |
                      (x & ((1 << boundary_bit) - 1));
        uint64_t xy = x & ((uint64_t) (-1) - ((1 << boundary_bit) - 1)) |
                      (y & ((1 << boundary_bit) - 1));
        int mutation_bit_xy = bit_distr(gen);
        int mutation_bit_yx = bit_distr(gen);
        xy ^= (1 << mutation_bit_xy);
        yx ^= (1 << mutation_bit_yx);
        return {xy, yx};
    }

    double genetic_opt(const std::function<double(double)> &f, double left, double right, double eps_f, int max_iters) {
        iters_amount = 0;

        const int selection_threshold = 8;
        std::mt19937 gen(42);
        std::uniform_real_distribution<double> population_distr(left, right);
        std::uniform_int_distribution<int> bit_distr(0, 63);
        std::vector<double> population(64);

        for (int i = 0; i < population.size(); i++) {
            population[i] = population_distr(gen);
        }

        double prev_best_value = f(left);
        for (int i = 0; i < max_iters; i++) {
            std::sort(population.begin(), population.end(),
                      [&f](double x, double y) { return f(x) < f(y); });

            double cur_best_value = f(population[0]);
            if (fabs(cur_best_value - prev_best_value) < eps_f) {
                break;
            }
            prev_best_value = cur_best_value;

            std::vector<double> new_population;
            for (int j = 0; j < selection_threshold; j++) {
                for (int j_pair = j + 1; j_pair < selection_threshold; j_pair++) {
                    uint64_t x = double_to_bitcode(population[j]);
                    uint64_t y = double_to_bitcode(population[j_pair]);
                    auto new_gen = cross_and_mutate(x, y, bit_distr, gen);
                    double xy = bitcode_to_double(new_gen.first);
                    double yx = bitcode_to_double(new_gen.second);

                    if (xy < left || xy > right || yx < left || yx > right) {
                        j_pair--;
                        continue;
                    }

                    new_population.push_back(xy);
                    new_population.push_back(yx);
                }
            }

            population = new_population;

            iters_amount++;
        }

        return *std::min_element(population.begin(), population.end(),
                                 [&f](double x, double y) { return f(x) < f(y); });
    }

}

double rastrigin_function(double x) {
    return 10 + x * x - 10 * cos(2 * M_PI * x);
}

double rastrigin_function_inv(double x) {
    return -(10 + x * x - 10 * cos(2 * M_PI * x));
}

void print_results(int iters_amount, double arg, double value, std::string method_name) {
    std::cout << "Метод " << method_name << ":" << std::endl;
    std::cout << "Число итераций: " << iters_amount << std::endl;
    std::cout << "Точка, к которой сошелся метод: " << arg << std::endl;
    std::cout << "Значение функции, к которому сошелся метод: " << value << std::endl;
}

int main() {
    const double left = -10, right = 30;
    const int max_iters = 10000;
    const double eps_a = 0.000001, eps_f = 0.000001;

    double ans = global_opt::random_search_opt(rastrigin_function, left, right, eps_f, max_iters);
    print_results(global_opt::iters_amount, ans, rastrigin_function(ans), "random search optimization");

    ans = global_opt::lipo_opt(rastrigin_function_inv, left, right, 123, eps_a, eps_f, max_iters);
    print_results(global_opt::iters_amount, ans, rastrigin_function(ans), "LIPO optimization");

    ans = global_opt::genetic_opt(rastrigin_function, left, right, eps_f, max_iters);
    print_results(global_opt::iters_amount, ans, rastrigin_function(ans), "genetic optimization");

    return 0;
}