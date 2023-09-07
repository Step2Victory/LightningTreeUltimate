#pragma once
#include <ranges>
#include <numeric>

namespace physics {

struct Vector {
    double x;
    double y;
    double z;
};

struct Vertex {
    double q;
    double Q;
    Vector r;
};

double distance(const Vector& lhs, const Vector& rhs);

double potential(const Vector& vector, const Vertex& vertex, double h);

template <typename T>
concept VertexSequence = std::ranges::range<T> && std::is_same_v<typename T::value_type, Vertex>;

template <VertexSequence Range>
double potential(const Vector& vector, const Range& range, double h) {

    return std::accumulate(std::ranges::begin(range), std::ranges::end(range), 0,
                           [h, vector](double ans, const Vertex& vertex) {
                               return ans + potential(vector, vertex, h);
                           });
}

}  // namespace physics