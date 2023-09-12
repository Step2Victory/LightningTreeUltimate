#include "constants.h"
#include "physics.h"
#include <cmath>

namespace physics {

double distance(const Vector& first, const Vector& second) {
    return std::sqrt((first.x - second.x) * (first.x - second.x) +
                     (first.y - second.y) * (first.y - second.y) +
                     (first.z - second.z) * (first.z - second.z));
}

double potential(const Vector& vector, const Vertex& vertex, double h) {
    static const double k = 1 / (4 * std::numbers::pi * EPSILON_0);
    double dist = distance(vertex.r, vector);

    if (dist < EPS) {
        return k * (vertex.q + vertex.Q) / h;
    }

    return k * (vertex.q + vertex.Q) / dist;
}

}  // namespace physics