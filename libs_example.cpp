#include <iostream>

#include "physics.h"

void physicsExample() {
    using namespace physics;
    Vertex vertex{
        .q = 0.0000001,
        .Q = 0.00001,
        .r =
            {
                .x = 0,
                .y = 100,
                .z = 10000,
            },
    };
    Vector r = {
        .x = 0,
        .y = 100,
        .z = 0,
    };
    auto ans = potential(r, vertex, 100);
    std::cout << "Example: " << ans << std::endl;
}

int main() {
    physicsExample();
}
