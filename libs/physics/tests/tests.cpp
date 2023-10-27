#include <gtest/gtest.h>
#include "physics.h"
#include <cmath>

using namespace physics;

TEST(PhysicsTest, potential_test_1) {
    physics::Vertex vertex{
        .q = 0.0000001,
        .Q = 0.00001,
        .r =
            {
                .x = 0,
                .y = 100,
                .z = 0,
            },
    };
    physics::Vector r = {
        .x = 0,
        .y = 100,
        .z = 0,
    };
    auto ans = physics::potential(r, vertex, 100);
    EXPECT_NEAR(ans, 907.74273101837866, 1e-7);
}

TEST(PhysicsTest, distance_test_1) {
    EXPECT_NEAR(distance({0, 0, 0}, {0, 0, 0}), 0, 1e-9);
    EXPECT_NEAR(distance({1, 0, 0}, {0, 0, 0}), 1, 1e-9);
    EXPECT_NEAR(distance({0, -2, 0}, {0, 2, 0}), 4, 1e-9);
    EXPECT_NEAR(distance({1, 0, 1}, {0, -1, 0}), std::sqrt(3), 1e-9);
}

// TODO: Add some other tests