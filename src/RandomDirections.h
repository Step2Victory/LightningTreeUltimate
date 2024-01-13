#pragma once
#include <random>
#include <vector>
#include <array>

std::vector<std::array<int, 3>> randomDirections() {
    static std::mt19937 gen(42);
    // static std::uniform_real_distribution<> dis(0, 26);
    // std::array<std::array<std::array<bool, 3>, 3>, 3> used;
    // for (int i = 0; i < 3; i++) {
    //     for (int j = 0; j < 3; j++) {
    //         for (int k = 0; k < 3; k++) {
    //             used[i][j][k] = false;
    //         }
    //     }
    // }
    // std::vector<std::array<int, 3>> directions;
    // int counter = 26;

    // while (counter) {
    //     int num = dis(gen);
    //     int i = num % 3;
    //     num /= 3;
    //     int j = num % 3;
    //     num /= 3;
    //     int k = num % 3;
    //     if (!used[i][j][k]) {
    //         directions.push_back({i, j, k});
    //         used[i][j][k] = true;
    //         counter--;
    //     }
    // }
    std::vector<std::array<int, 3>> directions;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                directions.push_back({i, j, k});
            }
        }
    }
    std::ranges::shuffle(directions, gen);
    return directions;
}


std::array<std::pair<double, std::array<int, 3>>, 26> sortDir(const std::array<std::pair<double, std::array<int, 3>>, 26>& potencial_dir){
    auto result_sort = potencial_dir;
    std::sort(result_sort.begin(), result_sort.end(), 
                [](const std::pair<double, std::array<int, 3>> &x,
                    const std::pair<double, std::array<int, 3>> &y) 
                    {
                        return std::abs(y.first) < std::abs(x.first) || (std::abs(y.first) < std::abs(x.first) && y.second[0] < x.second[0]);
                    });
    return result_sort;
}