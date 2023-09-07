#include <gtest/gtest.h>
#include "physics.h"

TEST(PhysicsTest, potential) {
    physics::Vertex vertex{
        .q = 0.0000001,
        .Q = 0.00001,
        .r =
            {
                .x = 0,
                .y = 100,
                .z = 10000,
            },
    };
    physics::Vector r = {
        .x = 0,
        .y = 100,
        .z = 0,
    };
    auto ans = physics::potential(r, vertex, 100);
    EXPECT_NEAR(ans, 0.1, 0.0000001);
}

TEST(PhysicsTest, distance) {
    // TO DO
}

// Add some other tests
// TO DO