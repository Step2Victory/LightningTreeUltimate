#include "LightningTree.h"
#include <numbers>
#include <random>
#include <cmath>

std::mt19937 gen(42);
std::uniform_real_distribution<> dis(0, 1);

constexpr double epsilon_0 = 8.854187817619999806e-12;
constexpr double kEps = 1e-9;

LightningTree::LightningTree(
    double h, double delta_t, double r, double R,
    size_t periphery_size, double q_plus_max, double q_minus_max,
    double Q_plus_s, double Q_minus_s, double resistance,
    double E_plus, double E_minus, double alpha, double beta,
    double sigma,
    std::function<double(const std::array<double, 3>&)>
        external_field_potential)
    : h(h),
      delta_t(delta_t),
      r(r),
      R(R),
      periphery_size(periphery_size),
      q_plus_max(q_plus_max),
      q_minus_max(q_minus_max),
      Q_plus_s(Q_plus_s),
      Q_minus_s(Q_minus_s),
      resistance(resistance),
      E_plus(E_plus),
      E_minus(E_minus),
      alpha(alpha),
      beta(beta),
      sigma(sigma),
      external_field_potential(external_field_potential) {
    // TO DO
}

void LightningTree::NextIter() {
    /*
    Описание метода
    */
    CountPotential();
    CountSigma();
    CountCurrent();
    double grow_time_period = 0.001;
    int grow_iter_period = grow_time_period / delta_t;
    if (iter_number % grow_iter_period == 0) {
        Grow();
        Delete();
    }
}

void LightningTree::CountPotential() {
    /*
     * Описание метода
     */
    for (auto& point : vertices) {
        double Phi = 0;
        for (auto& vertex : vertices) {
            double l = countDistance(vertex.coords, point.coords);
            double mirror_l = std::abs(std::sqrt(
                     std::pow((vertex.coords[0] - point.coords[0]), 2) +
                     std::pow((vertex.coords[1] - point.coords[1]), 2) +
                     std::pow((vertex.coords[2] + point.coords[2]), 2)));
            if (l < kEps) {
                Phi += vertex.q / (4 * std::numbers::pi * epsilon_0) *
                           (1 / (h + r) + 1 / (mirror_l + r)) +
                       vertex.Q / (4 * std::numbers::pi * epsilon_0) *
                           (1 / (h + R) + 1 / (mirror_l + R));
            }
            Phi += vertex.q / (4 * std::numbers::pi * epsilon_0) *
                       (1 / (l + r) + 1 / (mirror_l + r)) +
                   vertex.Q / (4 * std::numbers::pi * epsilon_0) *
                       (1 / (l + R) + 1 / (mirror_l + R));
        }
        point.Phi = Phi + external_field_potential(point.coords);
    }
}

double LightningTree::CountElectricity(const size_t v_from_id,
                                       const size_t v_to_id) const {
    /*
    Описание метода
    */
    double l = countDistance(vertices[v_from_id].coords, vertices[v_to_id].coords);
    double phi_from = vertices[v_from_id].Phi;
    double phi_to = vertices[v_to_id].Phi;

    // if (std::isinf((phi_from - phi_to) / l))
    // {
    //     throw std::runtime_error{ "Electric field along edge is
    //     infinity!" };
    // }
    return (phi_from - phi_to) / l;
}

void LightningTree::CountSigma() {
    /*
    Описание метода
    */
    for (auto edge : edges) {
        double E = CountElectricity(edge.from, edge.to);
        edge.sigma =
            edge.sigma * std::exp((alpha * E * E - beta) * delta_t);

        // if (std::isinf(sigma))
        // {
        //     throw std::runtime_error{ "Sigma is infinity!" };
        // }
    }
}

void LightningTree::CountCurrent() {
    /*
    Описание метода
    */
    for (auto edge : edges) {
        edge.current = std::numbers::pi * r * r * edge.sigma *
                       CountElectricity(edge.from, edge.to);
    }
}

void LightningTree::countCoords(std::array<double, 3>& result,
                                const size_t vertex_id,
                                const std::vector<int>& point) {
    result = {vertices[vertex_id].coords[0] + (1 - point[0]) * h,
              vertices[vertex_id].coords[1] + (1 - point[1]) * h,
              vertices[vertex_id].coords[2] + (1 - point[2]) * h};
}

double LightningTree::countDistance(const std::array<double, 3>& r_point, const std::array<double, 3>& l_point) const 
{
    return std::abs(std::sqrt(std::pow((r_point[0] - l_point[0]), 2) + std::pow((r_point[1] - l_point[1]), 2) + std::pow((r_point[2] - l_point[2]), 2)));
}

LightningTree::cubic_grid LightningTree::CreateNode(size_t edge_id, const std::vector<int>& point) {
    /*
    Описание метода
    */
    cubic_grid node;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                node[i][j][k] = -1;
            }
        }
    }
    node[1][1][1] = 0;
    vertices_peripherality.push_back(true);
    vertices_activity.push_back(true);
    node[2 - point[0]][2 - point[1]][2 - point[2]] = edge_id;
    return node;
}

void LightningTree::Transport() {
    /*
    Описание метода
    */
    for (auto& elem : graph) {
        size_t vertex_id = elem[1][1][1];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    if (elem[i][j][k] == -1 &&
                        (i == 1 && j == 1 && k == 1))
                        continue;

                    if (edges[elem[i][j][k]].from == vertex_id) {
                        vertices[vertex_id].q -=
                            edges[elem[i][j][k]].current * delta_t;
                    } else {
                        vertices[vertex_id].q +=
                            edges[elem[i][j][k]].current * delta_t;
                    }
                }
            }
        }
        if (vertices[vertex_id].q > q_plus_max) {
            vertices[vertex_id].Q +=
                vertices[vertex_id].q - q_plus_max;
            vertices[vertex_id].q = q_plus_max;
        } else if (vertices[vertex_id].q < q_minus_max) {
            vertices[vertex_id].Q +=
                vertices[vertex_id].q - q_minus_max;
            vertices[vertex_id].q = q_minus_max;
        }
    }
    iter_number++;
}

void LightningTree::Grow() {
    /*
    Описание метода
    */
    // TO DO
    std::vector<cubic_grid> temp_graph;
    for (size_t v_id = graph.size(); v_id-- >= 0;) {
        bool notGrow = true;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {

                    if (graph[v_id][i][j][k] != -1 &&
                        (i == 1 && j == 1 && k == 1))
                        continue;

                    std::array<double, 3> coords = {0, 0, 0};
                    countCoords(coords, v_id,
                                std::vector{i, j, k});
                    // Проверка на наличие вершины в векторе вершин по
                    // координатам.
                    Vertex new_vertex = {.q = 0,
                                         .Q = 0,
                                         .Phi = 0,
                                         .coords = coords,
                                         .number_edges = 0,
                                         .growless_iter_number = 0};
                    vertices.push_back(new_vertex);

                    if (GrowthCriterion(v_id,
                                        vertices.size() - 1)) {
                        notGrow = false;
                        vertices_peripherality[v_id] = false;
                        vertices[v_id].number_edges++;
                        // Проверка на наличие ребра в векторе рёбер
                        // по from и to вершинам. Занулять остальные
                        // переменные ребра???
                        Edge new_edge = {.from = v_id,
                                         .to = vertices.size() - 1,
                                         .current = 0,
                                         .sigma = 0};
                        edges.push_back(new_edge);
                        edges_activity.push_back(true);
                        graph[v_id][i][j][k] = edges.size() - 1;
                        temp_graph.push_back(CreateNode(
                            edges.size() - 1, std::vector{i, j, k}));
                    } else {
                        vertices.pop_back();
                    }
                }
            }
        }
        if (notGrow) 
        {
            vertices[v_id].growless_iter_number++;
        }
    }
    graph.insert(graph.end(), temp_graph.begin(), temp_graph.end());
}

void LightningTree::Delete() {
    /*
    Описание метода
    */
    for (size_t v_id = 0; v_id < vertices.size(); v_id++) {
        if (DeletionCriterion(v_id)) {
            vertices_activity[v_id] = false;
            vertices[v_id].growless_iter_number = 0;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    for (int k = 0; k < 3; k++) {
                        if (graph[v_id][i][j][k] == -1 &&
                            (i == 1 && j == 1 && k == 1))
                            continue;
                        size_t e_id = graph[v_id][i][j][k];
                        edges_activity[e_id] =
                            false;
                        if (vertices[edges[e_id]
                                         .from]
                                .number_edges == 1) {
                            vertices_peripherality
                                [edges[e_id]
                                     .from] = true;
                        }
                        vertices[edges[e_id]
                                         .from]
                                .number_edges--;
                    }
                }
            }
        }
    }
}

bool LightningTree::GrowthCriterion(const size_t v_from_id,
                                    const size_t v_to_id) const {
    /*
    Описание метода
    */
    double probability = dis(gen);
    double E = CountElectricity(v_from_id, v_to_id);
    if (E > E_plus) {
        return (1 - std::exp(-std::pow(((E - E_plus) / E_plus),
                                       2.5))) > probability;
    } else if (-E > E_minus) {
        return (1 - std::exp(-std::pow(((-E - E_minus) / E_minus),
                                       2.5))) > probability;
    }
    return false;
}
bool LightningTree::DeletionCriterion(size_t vertex_id) const {
    /*
    Описание метода
    */
    return (vertices_peripherality[vertex_id] &&
            vertices[vertex_id].growless_iter_number >
                periphery_size);
}