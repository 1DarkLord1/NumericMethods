#include <functional>
#include <cmath>
#include <cassert>
#include <iostream>


struct GaussLegendreIntegration {
    constexpr static size_t POINTS_NUMBER = 10;
    constexpr static double w[POINTS_NUMBER] = {
            0.295524224714752, 0.2955242247147529, 0.2692667193099963, 0.2692667193099963,
            0.2190863625159820, 0.2190863625159820, 0.1494513491505806, 0.1494513491505806,
            0.0666713443086881, 0.0666713443086881};
    constexpr static double x[POINTS_NUMBER] = {
            -0.1488743389816312, 0.1488743389816312, -0.4333953941292472, 0.4333953941292472,
            -0.6794095682990244, 0.6794095682990244, -0.8650633666889845, 0.8650633666889845,
            -0.9739065285171717, 0.9739065285171717};


    double integrate(double a, double b, const std::function<double(double)>& f) {
        double res = 0;
        double alpha = (b - a) / 2;
        double beta = (b + a) / 2;

        for (size_t i = 0; i < POINTS_NUMBER; i++) {
            res += w[i] * f(alpha * x[i] + beta);
        }

        return beta * res;
    }
};

// Task 2
double owen_function(double h, double a) {
    auto f = [h](double x) {
        return exp(-0.5 * h * h * (1 + x * x)) / (1 + x * x);
    };

    return GaussLegendreIntegration().integrate(0, a, f) / (2 * M_PI);
}

// Task 3
double trapezoid_integration(double a, double b, const std::function<double(double)>& f) {
    const int N = 1000;
    double step = (b - a) / N;
    double res = 0;
    for (int i = 1; i <= N; i++) {
        double x = a + (i - 1) * step;
        double y = a + i * step;
        res += 0.5 * (f(x) + f(y)) * (y - x);
    }

    return res;
}

double T(double x, double alpha, double t) {
    const double range = 100;
    auto f = [x, alpha, t](double y) {
        return exp(-x * x) * exp(-(x - y) * (x - y) / (2 * alpha * t));
    };

    return 0.5 * trapezoid_integration(-range, range, f) / (M_PI * alpha * t);
}

std::pair<double, double> find_min_max_T(double start_x, double end_x, double x_step, double alpha, double t) {
    double max_T = T(start_x, alpha, t), min_T = T(start_x, alpha, t);
    double x = start_x;

    while(x <= end_x) {
        max_T = std::max(max_T, T(x, alpha, t));
        min_T = std::min(min_T, T(x, alpha, t));
        x += x_step;
    }

    return {min_T, max_T};
}

double find_t(double start_x, double end_x) {
    double x_step = 0.01;
    double t_step = 0.5;
    double eps = 0.0001;
    double t = 0.5;

    while(true) {
        double min_T, max_T;
        std::tie(min_T, max_T) = find_min_max_T(start_x, end_x, x_step, 200, t);
        if (std::abs(max_T - min_T) <= eps) {
            break;
        }
        t += t_step;
    }

    return t;
}

int main() {
    double eps = 0.0001;

    double answers[] = {owen_function(1.34, 1), owen_function(0.21, 2), owen_function(2.55, 3.33)};
    double true_answers[] = {0.0410003, 0.169366, 0.00269307};
    for (size_t i = 0; i < 3; i++) {
        assert(std::abs(answers[i] - true_answers[i]) < eps);
    }
    std::cout << "Tests for Owen's function passed" << std::endl;

    std::cout << "t: " << find_t(1, 2) << std::endl;
    return 0;
}