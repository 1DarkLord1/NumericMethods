// Task 3

#include <functional>
#include <cmath>
#include <random>
#include <set>
#include <iostream>


namespace equation_solving {
    int iters_amount = 0;

    bool has_similar_root(const std::set<double> &roots, double root, double eps_a) {
        return std::any_of(roots.begin(), roots.end(),
                           [root, eps_a](double x) { return fabs(root - x) < eps_a; });
    }

    std::set<double>
    random_search_solving(const std::function<double(double)> &f, double left, double right, double eps_f,
                          double eps_sim, int roots_count, int max_iters) {
        iters_amount = 0;
        std::mt19937 gen(42);
        std::uniform_real_distribution<double> distr(left, right);
        std::set<double> roots;
        double x0 = (left + right) / 2;

        for (int i = 0; i < max_iters; i++) {
            if (roots.size() >= roots_count) {
                break;
            }

            if (fabs(f(x0)) < eps_f) {
                if (!has_similar_root(roots, x0, eps_sim)) {
                    roots.insert(x0);
                }
            }

            double delta = distr(gen);
            double pt = x0 + delta;

            if (pt < left || pt > right) {
                i--;
                continue;
            }

            x0 = pt;
            iters_amount++;
        }

        return roots;
    }

    double bisection(
            const std::function<double(double)> &f,
            double left,
            double right,
            double eps_f,
            double eps_a
    ) {
        while (fabs(right - left) >= eps_a && fabs(f(left)) >= eps_f) {
            auto mid = (right + left) / 2;
            auto val = f(mid);
            if (val <= 0) {
                left = mid;
            } else {
                right = mid;
            }

            iters_amount++;
        }

        return left;
    }

    std::set<double> hybrid_solving(
            const std::function<double(double)> &f,
            double left,
            double right,
            int segment_number,
            double eps_sim,
            double eps_f,
            double eps_a) {
        iters_amount = 0;
        double segment_len = (right - left) / segment_number;
        std::set<double> roots;

        for (int i = 0; i < segment_number; i++) {
            double cur_left = left + i * segment_len;
            double cur_right = left + (i + 1) * segment_len;
            double value_left = f(cur_left);
            double value_right = f(cur_right);

            if (value_left > 0 && value_right > 0 || value_left < 0 && value_right < 0) {
                continue;
            }

            double root_pt = bisection(f, cur_left, cur_right, eps_f, eps_a);
            if (!has_similar_root(roots, root_pt, eps_sim)) {
                roots.insert(root_pt);
            }
        }

        return roots;
    }
}

int main() {
    auto roots = equation_solving::random_search_solving(sin, -10, 10, 0.00001, 2, 7, 10000000);
    std::cout << "Random search method: " << std::endl;
    std::cout << "Iterations count: " << equation_solving::iters_amount << std::endl;
    for (auto root: roots) {
        std::cout << "x: " << root << ", " << "f(x): " << sin(root) << std::endl;
    }

    std::cout << std::endl;

    roots = equation_solving::hybrid_solving(sin, -10, 10, 10000, 2, 0.00001, 0.000001);
    std::cout << "Hybrid solving method: " << std::endl;
    std::cout << "Iterations count: " << equation_solving::iters_amount << std::endl;
    for (auto root: roots) {
        std::cout << "x: " << root << ", " << "f(x): " << sin(root) << std::endl;
    }
}