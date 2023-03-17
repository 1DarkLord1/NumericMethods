// Tasks 1 and 2

#include <cmath>
#include <cassert>
#include <iostream>
#include <random>
#include <functional>
#include <vector>

using vd = std::vector<double>;

namespace opt {
    int iters_counter = 0;

    double norm(const vd& v) {
        assert(v.size() == 2);
        return sqrt(v[0] * v[0] + v[1] * v[1]);
    }

    vd operator- (const vd& v1, const vd& v2) {
        assert(v1.size() == 2 && v2.size() == 2);
        return {v1[0] - v2[0], v1[1] - v2[1]};
    }

    vd operator+ (const vd& v1, const vd& v2) {
        assert(v1.size() == 2 && v2.size() == 2);
        return {v1[0] + v2[0], v1[1] + v2[1]};
    }

    vd operator/ (const vd& v1, double c) {
        assert(v1.size() == 2);
        return {v1[0] / c, v1[1] / c};
    }

    vd operator* (const std::vector<vd>& v1, const vd& v2) {
        assert(v1.size() == 2 && v2.size() == 2);
        return {v1[0][0] * v2[0] + v1[0][1] * v2[1], v1[1][0] * v2[0] + v1[1][1] * v2[1]};
    }

    vd bisection(
        const std::function<vd(vd)> &f,
        vd left,
        vd right,
        double epsilon_f,
        double epsilon_a
    ) {
        assert(left.size() == 2 && right.size() == 2);
        iters_counter = 0;
        while (norm(right - left) >= epsilon_a && norm(f(left)) >= epsilon_f) {
            auto mid = (right + left) / 2;
            auto val = f(mid);
            if (val[0] <= 0) {
                left[0] = mid[0];
            }
            else {
                right[0] = mid[0];
            }

            if (val[1] <= 0) {
                left[1] = mid[1];
            }
            else {
                right[1] = mid[1];
            }
            iters_counter++;
        }

        return left;
    }

    vd newton(
        const std::function<vd(vd)> &grad,
        const std::function<std::vector<vd>(vd)> &hessian_inv,
        vd x0,
        double epsilon_f,
        double epsilon_a
    ) {
        assert(x0.size() == 2);
        vd x = x0;
        iters_counter = 0;
        do {
            vd temp = x0;
            x0 = x;
            x = temp - hessian_inv(temp) * grad(temp);
            iters_counter++;
        } while(norm(x - x0) >= epsilon_a && norm(grad(x0)) >= epsilon_f);

        return x0;
    }
}

vd rastrigin_one_dim_grad(vd point) {
    // Игнорируем вторую координату
    return {2 * point[0] + 20 * M_PI * sin(2 * M_PI * point[0]), 0};
}

std::vector<vd> rastrigin_one_dim_hessian_inv(vd point) {
    // Игнорируем вторую координату
    return {{1.0 / (2 + 40 * M_PI * M_PI * cos(2 * M_PI * point[0])), 0}, {0, 0}};
}

vd rastrigin_two_dim_grad(vd point) {
    return {2 * point[0] + 20 * M_PI * sin(2 * M_PI * point[0]), 2 * point[1] + 20 * M_PI * sin(2 * M_PI * point[1])};
}

std::vector<vd> rastrigin_two_dim_hessian_inv(vd point) {
    double a = 2 + 40 * M_PI * M_PI * cos(2 * M_PI * point[0]);
    double d = 2 + 40 * M_PI * M_PI * cos(2 * M_PI * point[1]);
    return {{1.0 / a, 0}, {0, 1.0 / d}};
}

int main() {
    // Сгенерируем N случайных точек и посчитаем среднюю скорость сходимости.
    const int N = 1000;
    std::mt19937 gen(42);
    std::uniform_real_distribution<> dis(-10000, 10000);
    double eps_f = 0.0001, eps_a = 0.000001;

    // Сначала проверим реализацию на одномерном случае.
    // Метод бисекции
    double average_iters_amount = 0;
    int extremums_amount = 0;
    for (int i = 0; i < N; i++) {
        vd left = {dis(gen), 0}, right = {dis(gen), 0};
        if(left[0] > right[0]) {
            std::swap(left[0], right[0]);
        }
        vd ans = opt::bisection(rastrigin_one_dim_grad, left, right, eps_f, eps_a);
        if (fabs(rastrigin_one_dim_grad(ans)[0]) > eps_f) {
            // Метод сошелся к точке, значение в которой сильно отличается от 0
            continue;
        }
        average_iters_amount += opt::iters_counter;
        extremums_amount++;
    }

    average_iters_amount /= extremums_amount;
    std::cout << "(Одномерная функция Растригина, метод бисекции) Среднее число итераций: " << average_iters_amount <<
        " итераций, кол-во раз когда метод (почти) попал в экстремум: " << extremums_amount << " из " << N << std::endl;

    // Метод Ньютона
    average_iters_amount = 0;
    extremums_amount = 0;
    for (int i = 0; i < N; i++) {
        vd x0 = {dis(gen), dis(gen)};
        vd ans = opt::newton(rastrigin_one_dim_grad, rastrigin_one_dim_hessian_inv, x0, eps_f, eps_a);
        if (fabs(rastrigin_one_dim_grad(ans)[0]) > eps_f) {
            // Метод сошелся к точке, значение в которой сильно отличается от 0
            continue;
        }
        average_iters_amount += opt::iters_counter;
        extremums_amount++;
    }

    average_iters_amount /= extremums_amount;
    std::cout << "(Одномерная функция Растригина, метод Ньютона) Среднее число итераций: " << average_iters_amount <<
              " итераций, кол-во раз когда метод (почти) попал в экстремум: " << extremums_amount << " из " << N << std::endl;

    // Теперь проверим реализацию на двумерном случае.
    // Метод бисекции
    average_iters_amount = 0;
    extremums_amount = 0;
    for (int i = 0; i < N; i++) {
        vd left = {dis(gen), dis(gen)}, right = {dis(gen), dis(gen)};
        if(left[0] > right[0]) {
            std::swap(left[0], right[0]);
        }
        if(left[1] > right[1]) {
            std::swap(left[1], right[1]);
        }
        vd ans = opt::bisection(rastrigin_two_dim_grad, left, right, eps_f, eps_a);
        if (opt::norm(rastrigin_two_dim_grad(ans)) > eps_f) {
            // Метод сошелся к точке, значение в которой сильно отличается от 0
            continue;
        }
        average_iters_amount += opt::iters_counter;
        extremums_amount++;
    }

    average_iters_amount /= extremums_amount;
    std::cout << "(Двумерная функция Растригина, метод бисекции) Среднее число итераций: " << average_iters_amount <<
              " итераций, кол-во раз когда метод (почти) попал в экстремум: " << extremums_amount << " из " << N << std::endl;

    // Метод Ньютона
    average_iters_amount = 0;
    extremums_amount = 0;
    for (int i = 0; i < N; i++) {
        vd x0 = {dis(gen), dis(gen)};
        vd ans = opt::newton(rastrigin_two_dim_grad, rastrigin_two_dim_hessian_inv, x0, eps_f, eps_a);
        if (opt::norm(rastrigin_two_dim_grad(ans)) > eps_f) {
            // Метод сошелся к точке, значение в которой сильно отличается от 0
            continue;
        }
        average_iters_amount += opt::iters_counter;
        extremums_amount++;
    }

    average_iters_amount /= extremums_amount;
    std::cout << "(Двумерная функция Растригина, метод Ньютона) Среднее число итераций: " << average_iters_amount <<
              " итераций, кол-во раз когда метод (почти) попал в экстремум: " << extremums_amount << " из " << N << std::endl;

    // Метод бисекции сходится примерно за 34 итерации в обоих случаях.
    // При допустимых погрешностях eps_f = 0.0001, eps_a = 0.000001 сходится к более менее хорошей точке довольно нечасто
    // (224 раза из 1000 и 96 раз из 1000).
    // С методом Ньютона возникли проблемы, почти всегда он сходится за 2 итерации и сходится
    // к очень плохой точке.

    return 0;
}
