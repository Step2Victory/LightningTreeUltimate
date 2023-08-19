#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <numbers>

constexpr double epsilon_0 = 8.854187817619999806e-12;
constexpr double kEps = 1e-9;

struct ChargeLayer {
    double p_0;
    double h;
    double L;
    std::array<double, 3> r;
    double a;
};

double Potential(const std::array<double, 3>& point,
                 const std::vector<std::vector<std::vector<double>>>& q_values,
                 const std::array<double, 3>& start, double h) {
    double result = 0;
    for (int i = 0; i < q_values.size(); i++) {
        for (int j = 0; j < q_values[0].size(); j++) {
            for (int k = 0; k < q_values[0][0].size(); k++) {
                std::array<double, 3> r = {start[0] + h * i, start[1] + h * j,
                                           start[2] + h * k};

                double l = std::abs(std::sqrt(std::pow((r[0] - point[0]), 2) +
                                              std::pow((r[1] - point[1]), 2) +
                                              std::pow((r[2] - point[2]), 2)));
                double mirror_l =
                    std::abs(std::sqrt(std::pow((r[0] - point[0]), 2) +
                                       std::pow((r[1] - point[1]), 2) +
                                       std::pow((r[2] + point[2]), 2)));

                if (l < kEps) {
                    result += q_values[i][j][k] /
                              (4 * std::numbers::pi * epsilon_0) *
                              (1 / (h)-1 / (mirror_l));
                }
                result += q_values[i][j][k] /
                          (4 * std::numbers::pi * epsilon_0) *
                          (1 / (l)-1 / (mirror_l));
            }
        }
    }
    return result;
}

std::function<double(const std::array<double, 3>&)> constExternalField(
    double external_field) {
    std::function<double(const std::array<double, 3>&)> result =
        [external_field =
             std::move(external_field)](const std::array<double, 3>& coords) {
            return coords[2] * external_field;
        };
    return result;
}

std::function<double(const std::array<double, 3>&)> countExternalField(
    const std::vector<ChargeLayer>& layers, const std::array<double, 3>& start,
    const std::array<double, 3>& end, double h) {
    std::vector<std::vector<std::vector<double>>> q_values(
        (end[0] - start[0]) / h,
        std::vector<std::vector<double>>(
            (end[1] - start[1]) / h,
            std::vector<double>((end[2] - start[2]) / h)));

    std::vector<std::vector<std::vector<double>>> potential_values(q_values);

    for (int i = 0; i < q_values.size(); i++) {
        for (int j = 0; j < q_values[0].size(); j++) {
            for (int k = 0; k < q_values[0][0].size(); k++) {
                double q = 0;
                std::array<double, 3> point = {
                    start[0] + h * i, start[1] + h * j, start[2] + h * k};
                for (auto& layer : layers) {
                    std::array<double, 3> cur = {point[0] - layer.r[0],
                                                 point[1] - layer.r[1],
                                                 point[2] - layer.r[2]};

                    auto expr =
                        -std::pow(std::sqrt(cur[2] * cur[2]) / layer.h,
                                  2 * layer.a) -
                        std::pow(std::sqrt(cur[0] * cur[0] + cur[1] * cur[1]) /
                                     layer.L,
                                 2 * layer.a);
                    q += layer.p_0 * std::exp(expr) * h * h * h;
                }
                q_values[i][j][k] = q;
            }
        }
    }
    for (int i = 0; i < q_values.size(); i++) {
        for (int j = 0; j < q_values[0].size(); j++) {
            for (int k = 0; k < q_values[0][0].size(); k++) {
                std::array<double, 3> point = {
                    start[0] + h * i, start[1] + h * j, start[2] + h * k};
                potential_values[i][j][k] =
                    Potential(point, q_values, start, h);
            }
        }
    }

    std::function<double(const std::array<double, 3>&)> result =
        [potential_values = std::move(potential_values), r = std::move(start),
         h = std::move(h)](const std::array<double, 3>& coords) {
            std::array<int, 3> shift = {
                static_cast<int>((coords[0] - r[0]) / h),
                static_cast<int>((coords[1] - r[1]) / h),
                static_cast<int>((coords[2] - r[2]) / h)};
            // std::array<int, 3> shift = {static_cast((coords[0] - r[0])/h),
            // static_cast((coords[1] - r[1])/h), static_cast((coords[2] -
            // r[2])/h)};
            if (shift[0] < 0 || shift[1] < 0 || shift[1] < 0 ||
                shift[0] >= potential_values.size() ||
                shift[1] >= potential_values[0].size() ||
                shift[1] >= potential_values[0][0].size()) {
                std::cout << "Выход за границу расчетной области!\n";
                // throw std::runtime_error{"Выход за границу расчетной
                // области!"};
            }
            return potential_values[shift[0]][shift[1]][shift[2]];
        };
    return result;
}