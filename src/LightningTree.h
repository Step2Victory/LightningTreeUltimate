#include <vector>
#include <array>
#include <functional>
#include <filesystem>
#include <random>
#include <unordered_map>
#include <map>


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
                  double sigma, std::array<double, 3> start_r, std::array<double, 3> end_r,
                  double degree_probability_growth, int seed);
    void NextIter();
    void CountSigma();
    double Potential(const std::array<double, 3>&);
    void CountPotential();
    double CountElectricity(size_t, size_t) const;
    void CountCurrent();
    std::array<double, 3> countCoords(size_t, const std::array<int, 3>&);
    std::array<int, 3> countInternalCoords(size_t, const std::array<int, 3>&);

    double countDistance(const std::array<double, 3>&, const std::array<double, 3>&) const;
    cubic_grid CreateNode(size_t, const std::array<int, 3>&);

    void Transport();
    void Grow();
    void Delete();

    bool GrowthCriterion(size_t, size_t) const;
    bool DeletionCriterion(size_t) const;

    void AllParams() const;
    void Info() const;
    void WriteResponse(int response) const;
    void ReturnFiles(const std::filesystem::path&);
    void ReturnPhi(const std::filesystem::path&, const std::array<double, 3>&, const std::array<double, 3>&);

    const std::array<double, 3> start_r;
    const std::array<double, 3> end_r;

private:

    bool TryAddEdge(size_t v_from_id, const std::array<int, 3>& dir);

    struct Vertex {
        double q;
        double Q;
        double Phi;
        std::array<double, 3> coords;
        std::array<int, 3> internal_coords;
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

    struct Less {
        bool operator()(const std::array<int, 3>& lhs, const std::array<int, 3>& rhs) const {
            if (lhs[0] == rhs[0] && lhs[1] == rhs[1]) {
                return lhs[2] < rhs[2];
            }
            if (lhs[0] == rhs[0]) {
                return lhs[1] < rhs[1];
            }
            return lhs[0] < rhs[0];
        }
    };

    std::map<std::array<int, 3>, int, Less> internal_coords_to_id; 
    size_t iter_number;

    int seed;
    mutable std::mt19937 gen;
    mutable std::uniform_real_distribution<> dis;
};
