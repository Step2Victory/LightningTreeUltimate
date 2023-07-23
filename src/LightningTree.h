#include <vector>
#include <array>
#include <functional>
#include <filesystem>
#include <random>

class LightningTree {
private:
    using cubic_grid =
        std::array<std::array<std::array<int, 3>, 3>, 3>;

public:
    LightningTree(const std::filesystem::path& path_to_config_file);
    LightningTree(double h, double delta_t, double r, double R,
                  size_t periphery_size, double q_plus_max,
                  double q_minus_max, double Q_plus_s,
                  double Q_minus_s, double resistance, double E_plus,
                  double E_minus, double alpha, double beta,
                  double sigma,
                  std::function<double(const std::array<double, 3>&)>
                      external_field_potential, int seed);
    void NextIter();
    void CountSigma();
    void CountPotential();
    double CountElectricity(size_t, size_t) const;
    void CountCurrent();
    std::array<double, 3> countCoords(size_t, const std::array<int, 3>&);
    double countDistance(const std::array<double, 3>&, const std::array<double, 3>&) const;
    cubic_grid CreateNode(size_t, const std::array<int, 3>&);

    void Transport();
    void Grow();
    void Delete();

    bool GrowthCriterion(size_t, size_t) const;
    bool DeletionCriterion(size_t) const;

    void AllParams();
    void Info();
    void ReturnFiles(const std::filesystem::path&);
    void ReturnPhi(double, double);
private:
    struct Vertex {
        double q;
        double Q;
        double Phi;
        std::array<double, 3> coords;
        size_t number_edges;
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
    double degree_probability_growth;
    size_t periphery_size;
    std::function<double(const std::array<double, 3>&)>
        external_field_potential;

    std::vector<Vertex> vertices;
    std::vector<bool> vertices_peripherality;
    std::vector<bool> vertices_activity;
    std::vector<Edge> edges;
    std::vector<bool> edges_activity;

    //здесь индексы в массиве ребер

    std::vector<cubic_grid> graph;

    size_t iter_number;

    int seed;
    mutable std::mt19937 gen;
    mutable std::uniform_real_distribution<> dis;
};
