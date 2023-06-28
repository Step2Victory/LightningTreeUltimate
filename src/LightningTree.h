#include <vector>
#include <array>
#include <functional>

class LightningTree {

public:
    LightningTree(double h, double delta_t, double r, double R,
                  size_t periphery_size, double q_plus_max,
                  double q_minus_max, double Q_plus_s,
                  double Q_minus_s, double resistance, double E_plus,
                  double E_minus, double alpha, double beta,
                  double sigma,
                  std::function<double(double, double, double)>
                      external_field_potential);
    void NextIter();
    void CountSigma();
    void CountElectricity();
    void CountCurrent();

    void Transport();
    void Grow();
    void Delete();

    bool GrowthCriterion(size_t vertex_id,
                         const std::array<double, 3>& coords) const;
    bool DeletionCriterion(size_t vertex_id) const;

private:
    struct Vertex {
        double q;
        double Q;
        double E;
        std::array<double, 3> coords;
        size_t growless_iter_number;
    };

    struct Edge {
        size_t from;
        size_t to;
        double current;
        double sigma;
    };

    double h;
    double delta_t;
    double r;
    double R;
    double q_plus_max;
    double q_minus_max;
    double Q_plus_s;
    double Q_minus_s;
    double resistance;
    double E_plus;
    double E_minus;
    double alpha;
    double beta;
    double sigma;
    size_t periphery_size;
    std::function<double(double, double, double)>
        external_field_potential;

    std::vector<Vertex> vertices;
    std::vector<Edge> edges;

    //здесь индексы в массиве ребер
    using cubic_grid =
        std::array<std::array<std::array<int, 3>, 3>, 3>;
    std::vector<cubic_grid> graph;

    size_t iter_number;
};
