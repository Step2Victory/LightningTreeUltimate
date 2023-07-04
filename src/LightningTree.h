#include <vector>
#include <array>
#include <functional>
#include <cmath>

// ???
std::mt19937 gen(42);
std::uniform_real_distribution<> dis(0, 1);

constexpr double epsilon_0 = 8.854187817619999806e-12;
constexpr double kEps = 1e-9;
constexpr double PI = 3.14159265358979323846;

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
    void CountPotential();
    double CountElectricity(const size_t, const size_t) const;
    void CountCurrent();
    void countCoords(std::array<double, 3>&, const size_t, const std::vector<int>);
    cubic_grid CreateNode(size_t, size_t, const std::vector<int>&); // ???

    void Transport();
    void Grow();
    void Delete();

    bool GrowthCriterion(const size_t, const size_t) const;
    bool DeletionCriterion(size_t) const;

private:
    struct Vertex {
        double q;
        double Q;
        double Phi;
        std::array<double, 3> coords;
        size_t growless_iter_number; // ???
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
    double resistance; // ???
    double E_plus;
    double E_minus;
    double alpha;
    double beta;
    double sigma;
    size_t periphery_size; // ???
    std::function<double(double, double, double)>
        external_field_potential; // ???

    std::vector<Vertex> vertices;
    std::vector<Edge> edges;

    using cubic_grid =
        std::array<std::array<std::array<int, 3>, 3>, 3>;

    std::vector<cubic_grid> graph;
    
    size_t iter_number;
};