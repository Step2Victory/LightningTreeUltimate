#include <random>
#include <vector>
#include <array>

std::mt19937 gen(42);
std::uniform_real_distribution<> dis(0, 26);

std::vector<std::array<int, 3>> randomDirections()
{
    std::array<std::array<std::array<bool, 3>, 3>, 3> used = {{{false}}};
    std::vector<std::array<int, 3>> directions;
    int counter = 26;

    // for (size_t i = 0; i < 3; ++i)
    // {
    //     for (size_t j = 0; j < 3; ++j)
    //     {
    //         for (size_t k = 0; k < 3; ++k)
    //         {
    //             used[i][j][k] = false;
    //         }
    //     }
    // }

    while (counter)
    {
        int num = dis(gen);
        int i = num % 3;
        num /= 3;
        int j = num % 3;
        num /= 3;
        int k = num % 3;
        if (!used[i][j][k]) 
        {
            directions.push_back({i, j, k});
            used[i][j][k] = true;
            counter--;
        }
    }
    return directions;
}