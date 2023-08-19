#include "LightningTree.h"

#include "yaml-cpp/yaml.h"
#include <glog/logging.h>

#include <numbers>
#include <cmath>
#include <fstream>
#include <string>
#include <exception>

#include "ExternalField.h"
#include "RandomDirections.h"

LightningTree::LightningTree(const std::filesystem::path& path_to_config_file) {
    YAML::Node config = YAML::LoadFile(path_to_config_file);
    h = config["h"].as<double>();
    delta_t = config["delta_t"].as<double>();
    r = config["r"].as<double>();
    R = config["R"].as<double>();
    periphery_size = config["periphery_size"].as<double>();
    q_plus_max = config["q_plus_max"].as<double>();
    q_minus_max = config["q_minus_max"].as<double>();
    Q_plus_s = config["Q_plus_s"].as<double>();
    Q_minus_s = config["Q_minus_s"].as<double>();
    resistance = config["resistance"].as<double>();
    E_plus = config["E_plus"].as<double>();
    E_minus = config["E_minus"].as<double>();
    alpha = config["alpha"].as<double>();
    beta = config["beta"].as<double>();
    sigma = config["sigma"].as<double>();
    start_r = {config["start_x"].as<double>(), config["start_y"].as<double>(),
               config["start_z"].as<double>()};
    end_r = {
        config["end_x"].as<double>(),
        config["end_y"].as<double>(),
        config["end_z"].as<double>(),
    };
    std::vector<ChargeLayer> layers;
    for (const auto& layer : config["external_field_layers"]) {
        layers.push_back(ChargeLayer{.p_0 = layer["p_0"].as<double>(),
                                     .h = layer["h"].as<double>(),
                                     .L = layer["L"].as<double>(),
                                     .r = layer["r"].as<std::array<double, 3>>(),
                                     .a = layer["a"].as<double>()});
    }
    external_field_potential = countExternalField(layers, start_r, end_r, h);

    seed = config["seed"].as<int>();
    degree_probability_growth = config["degree_probability_growth"].as<double>();
    gen = std::mt19937(seed);
    dis = std::uniform_real_distribution<>(0, 1);
    iter_number = 0;
    vertices.push_back(Vertex{.q = 0,
                              .Q = 0,
                              .Phi = 0,
                              .coords = {(end_r[0] + start_r[0]) / 2, (end_r[1] + start_r[1]) / 2,
                                         (end_r[2] + start_r[2]) / 2 - h},
                              .internal_coords = {0, 0, 0},
                              .number_edges = 0,
                              .growless_iter_number = 0});
    vertices.push_back(Vertex{.q = 0,
                              .Q = 0,
                              .Phi = 0,
                              .coords = {(end_r[0] + start_r[0]) / 2, (end_r[1] + start_r[1]) / 2,
                                         (end_r[2] + start_r[2]) / 2},
                              .internal_coords = {0, 0, 1},
                              .number_edges = 0,
                              .growless_iter_number = 0});
    edges.push_back(Edge{.from = 0, .to = 1, .current = 0, .sigma = sigma});
    edges_activity.push_back(true);
    graph.push_back(CreateNode(edges.size() - 1, {1, 1, 2}));
    graph.push_back(CreateNode(edges.size() - 1, {1, 1, 0}));
}

LightningTree::LightningTree(double h, double delta_t, double r, double R, size_t periphery_size,
                             double q_plus_max, double q_minus_max, double Q_plus_s,
                             double Q_minus_s, double resistance, double E_plus, double E_minus,
                             double alpha, double beta, double sigma, std::array<double, 3> start_r,
                             std::array<double, 3> end_r, double degree_probability_growth,
                             int seed)
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
      start_r(start_r),
      end_r(end_r),
      //   external_field_potential(countExternalField(
      //     std::vector({
      //         ChargeLayer{.p_0 = -0.0000007, .h = h, .L = 200, .r =
      //         start_r, .a = 1}, ChargeLayer{.p_0 = 0.0000007, .h =
      //         h, .L = 200, .r = end_r, .a = 1}}
      //         ),
      //     std::array<double, 3>{end_r[0] - start_r[0] - 5 * h,
      //     end_r[1] - start_r[1] - 5 * h, end_r[2] - start_r[2] - 20
      //     * h}, std::array<double, 3>{end_r[0] - start_r[0] + 5 *
      //     h, end_r[1] - start_r[1] + 5 * h, end_r[2] - start_r[2] +
      //     20 * h}, h)
      //   ),
      external_field_potential(constExternalField(100000.0)),
      degree_probability_growth(degree_probability_growth),
      seed(seed) {
    gen = std::mt19937(seed);
    dis = std::uniform_real_distribution<>(0, 1);
    iter_number = 0;
    vertices.push_back(Vertex{.q = 0,
                              .Q = 0,
                              .Phi = 0,
                              .coords = {(end_r[0] + start_r[0]) / 2, (end_r[1] + start_r[1]) / 2,
                                         (end_r[2] + start_r[2]) / 2 - h},
                              .internal_coords = {0, 0, 0},
                              .number_edges = 0,
                              .growless_iter_number = 0});
    vertices.push_back(Vertex{.q = 0,
                              .Q = 0,
                              .Phi = 0,
                              .coords = {(end_r[0] + start_r[0]) / 2, (end_r[1] + start_r[1]) / 2,
                                         (end_r[2] + start_r[2]) / 2},
                              .internal_coords = {0, 0, 1},
                              .number_edges = 0,
                              .growless_iter_number = 0});
    edges.push_back(Edge{.from = 0, .to = 1, .current = 0, .sigma = sigma});
    edges_activity.push_back(true);
    graph.push_back(CreateNode(edges.size() - 1, {1, 1, 2}));
    graph.push_back(CreateNode(edges.size() - 1, {1, 1, 0}));
}

void LightningTree::NextIter() {
    CountPotential();
    CountSigma();
    CountCurrent();
    Transport();

    double grow_time_period = 0.001;
    int grow_iter_period = grow_time_period / delta_t;
    if (iter_number % grow_iter_period == 0) {
        Grow();
        Delete();
    }
}

double LightningTree::Potential(const std::array<double, 3>& coords) {
    double Phi = 0;
    for (auto& vertex : vertices) {
        double l = countDistance(vertex.coords, coords);
        double mirror_l = countDistance(vertex.coords, {coords[0], coords[1], -coords[2]});
        if (l < kEps) {
            Phi +=
                vertex.q / (4 * std::numbers::pi * epsilon_0) * (1 / (h + r) - 1 / (mirror_l + r)) +
                vertex.Q / (4 * std::numbers::pi * epsilon_0) * (1 / (h + R) - 1 / (mirror_l + R));
        } else {
            Phi +=
                vertex.q / (4 * std::numbers::pi * epsilon_0) * (1 / (l + r) - 1 / (mirror_l + r)) +
                vertex.Q / (4 * std::numbers::pi * epsilon_0) * (1 / (l + R) - 1 / (mirror_l + R));
        }
    }
    return Phi;
}

void LightningTree::CountPotential() {
    for (auto& point : vertices) {
        point.Phi = Potential(point.coords) + external_field_potential(point.coords);
    }
}

double LightningTree::CountElectricity(const size_t v_from_id, const size_t v_to_id) const {
    double l = countDistance(vertices[v_from_id].coords, vertices[v_to_id].coords);
    double phi_from = vertices[v_from_id].Phi;
    double phi_to = vertices[v_to_id].Phi;

    if (std::isinf((phi_from - phi_to) / l)) {
        LOG(INFO) << "Electric field in " << v_from_id << ": " << phi_from << ". "
                  << "Electric field in " << v_to_id << ": " << phi_to;
        throw std::runtime_error{"Electric field along edge is infinity !"};
    }

    return (phi_from - phi_to) / l;
}

void LightningTree::CountSigma() {
    for (auto& edge : edges) {
        double E = CountElectricity(edge.from, edge.to);
        edge.sigma = edge.sigma * std::exp((alpha * E * E - beta) * delta_t);

        if (std::isinf(sigma)) {
            throw std::runtime_error{"Sigma in edge is infinity!"};
        }
    }
}

void LightningTree::CountCurrent() {
    for (auto& edge : edges) {
        edge.current = std::numbers::pi * r * r * edge.sigma * CountElectricity(edge.from, edge.to);
    }
}

std::array<double, 3> LightningTree::countCoords(const size_t vertex_id,
                                                 const std::array<int, 3>& dir) {
    return {vertices[vertex_id].coords[0] + (dir[0] - 1) * h,
            vertices[vertex_id].coords[1] + (dir[1] - 1) * h,
            vertices[vertex_id].coords[2] + (dir[2] - 1) * h};
}

std::array<int, 3> LightningTree::countInternalCoords(size_t vertex_id,
                                                      const std::array<int, 3>& dir) {
    return {vertices[vertex_id].internal_coords[0] + (dir[0] - 1),
            vertices[vertex_id].internal_coords[1] + (dir[1] - 1),
            vertices[vertex_id].internal_coords[2] + (dir[2] - 1)};
}

double LightningTree::countDistance(const std::array<double, 3>& r_point,
                                    const std::array<double, 3>& l_point) const {
    return std::abs(std::sqrt(std::pow((r_point[0] - l_point[0]), 2) +
                              std::pow((r_point[1] - l_point[1]), 2) +
                              std::pow((r_point[2] - l_point[2]), 2)));
}

LightningTree::cubic_grid LightningTree::CreateNode(size_t edge_id, const std::array<int, 3>& dir) {

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
    node[2 - dir[0]][2 - dir[1]][2 - dir[2]] = edge_id;
    return node;
}

void LightningTree::Transport() {
    for (size_t g_id = 0; g_id < graph.size(); g_id++) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    int e_id = graph[g_id][i][j][k];
                    if (e_id == -1 || !edges_activity[e_id] || (i == 1 && j == 1 && k == 1))
                        continue;

                    if (g_id == edges[e_id].from) {
                        vertices[g_id].q -= edges[e_id].current * delta_t;
                    } else {
                        vertices[g_id].q += edges[e_id].current * delta_t;
                    }
                }
            }
        }
        if (vertices[g_id].q > q_plus_max) {
            vertices[g_id].Q += vertices[g_id].q - q_plus_max;
            vertices[g_id].q = q_plus_max;
        } else if (vertices[g_id].q < q_minus_max) {
            vertices[g_id].Q += vertices[g_id].q - q_minus_max;
            vertices[g_id].q = q_minus_max;
        }
    }
    iter_number++;
}

bool LightningTree::TryAddEdge(size_t v_from_id, const std::array<int, 3>& dir) {

    std::array<double, 3> coords = countCoords(v_from_id, {dir[0], dir[1], dir[2]});
    std::array<int, 3> internal_coords = countInternalCoords(v_from_id, {dir[0], dir[1], dir[2]});

    // Проверка на наличие вершины в векторе
    // вершин по координатам.
    size_t v_to_id;
    bool new_vertex = true;
    if (internal_coords_to_id.contains(internal_coords)) {
        v_to_id = internal_coords_to_id[internal_coords];
        new_vertex = false;
    }
    if (new_vertex) {
        Vertex new_vertex = {.q = 0,
                             .Q = 0,
                             .Phi = Potential(coords) + external_field_potential(coords),
                             .coords = coords,
                             .internal_coords = internal_coords,
                             .number_edges = 0,
                             .growless_iter_number = 0};
        vertices.push_back(new_vertex);
        v_to_id = vertices.size() - 1;
        internal_coords_to_id[internal_coords] = v_to_id;
    }

    if (GrowthCriterion(v_from_id, v_to_id)) {
        // std::cout << "grow" << '\n';
        vertices_peripherality[v_from_id] = false;
        vertices[v_from_id].number_edges++;

        Edge new_edge = {.from = v_from_id, .to = v_to_id, .current = 0, .sigma = sigma};
        edges.push_back(new_edge);
        edges_activity.push_back(true);
        graph[v_from_id][dir[0]][dir[1]][dir[2]] = edges.size() - 1;
        if (new_vertex) {
            graph.push_back(CreateNode(edges.size() - 1, dir));
            vertices_peripherality.push_back(true);
        } else {
            vertices_activity[v_to_id] = true;
            graph[v_to_id][2 - dir[0]][2 - dir[1]][2 - dir[2]] = edges.size() - 1;
        }
        return true;
    } else if (new_vertex) {
        vertices.pop_back();
    }
    return false;
}

void LightningTree::Grow() {
    std::vector<cubic_grid> temp_graph;
    int number_of_vertices_before_grow = graph.size();
    for (size_t v_from_id = number_of_vertices_before_grow; v_from_id-- > 0;) {
        if (!vertices_activity[v_from_id])
            continue;
        bool notGrow = true;
        auto directions = randomDirections();
        for (auto dir : directions) {
            int i = dir[0];
            int j = dir[1];
            int k = dir[2];
            int e_id = graph[v_from_id][i][j][k];

            if (i == 1 && j == 1 && k == 1)
                continue;

            if (e_id == -1) {
                notGrow = !TryAddEdge(v_from_id, dir);

            } else if (!edges_activity[e_id]) {
                // проверка неактивного ребра
                if (GrowthCriterion(v_from_id, edges[e_id].to)) {
                    // std::cout << "grow" << '\n';
                    notGrow = false;
                    vertices_peripherality[v_from_id] = false;
                    vertices[v_from_id].number_edges++;

                    vertices_activity[edges[e_id].to] = true;
                    edges_activity[e_id] = true;
                }

            } else
                continue;
        }
        if (notGrow) {
            vertices[v_from_id].growless_iter_number++;
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
            vertices[v_id].number_edges = 0;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    for (int k = 0; k < 3; k++) {
                        int e_id = graph[v_id][i][j][k];
                        if ((e_id == -1 || !edges_activity[e_id]) || (i == 1 && j == 1 && k == 1))
                            continue;
                        edges_activity[e_id] = false;
                        vertices[edges[e_id].from].number_edges--;

                        if (vertices[edges[e_id].from].number_edges == 1) {
                            vertices_peripherality[edges[e_id].from] = true;
                        }
                        if (vertices[edges[e_id].from].number_edges == 0) {
                            vertices_activity[edges[e_id].from] = true;
                        }
                    }
                }
            }
        }
    }
}

bool LightningTree::GrowthCriterion(const size_t v_from_id, const size_t v_to_id) const {
    /*
    Описание метода
    */
    // std::cout << "try grow" << '\n';
    double probability = dis(gen);
    double E = CountElectricity(v_from_id, v_to_id);
    if (E > E_plus) {
        return (1 - std::exp(-std::pow(((E - E_plus) / E_plus), degree_probability_growth))) >
               probability;
    } else if (-E > E_minus) {
        return (1 - std::exp(-std::pow(((-E - E_minus) / E_minus), degree_probability_growth))) >
               probability;
    }
    return false;
}

bool LightningTree::DeletionCriterion(size_t vertex_id) const {
    /*
    Описание метода
    */
    return (vertices_peripherality[vertex_id] &&
            vertices[vertex_id].growless_iter_number > periphery_size);
}

void LightningTree::AllParams() const {
    LOG(INFO) << "h: " << h << '\n'
              << "delta_t: " << delta_t << '\n'
              << "r: " << r << '\n'
              << "R: " << R << '\n'
              << "periphery_size: " << periphery_size << '\n'
              << "q_plus_max: " << q_plus_max << '\n'
              << "q_minus_max: " << q_minus_max << '\n'
              << "Q_plus_s: " << Q_plus_s << '\n'
              << "Q_minus_s: " << Q_minus_s << '\n'
              << "resistance: " << resistance << '\n'
              << "E_plus: " << E_plus << '\n'
              << "E_minus: " << E_minus << '\n'
              << "alpha: " << alpha << '\n'
              << "beta: " << beta << '\n'
              << "sigma: " << sigma << '\n'
              << '\n'
              << "seed: " << seed << '\n';
}

void LightningTree::Info() const {
    LOG(INFO) << "iter_number: " << iter_number << '\t' << "number_of_vertices: " << vertices.size()
              << '\t';
}

void LightningTree::WriteResponse(int response) const {
    std::cout << response << ' ' << iter_number << ' ' << vertices.size() << ' '
              << iter_number * delta_t << std::endl;
}

void LightningTree::ReturnFiles(const std::filesystem::path& path_to_data) {
    std::ofstream fout(path_to_data / "vertex_table.txt");
    fout << "id" << ' ' << 'q' << ' ' << 'Q' << ' ' << 'x' << ' ' << 'y' << ' ' << 'z' << ' '
         << "phi" << '\n';
    for (int i = 0; i < vertices.size(); i++) {
        fout << "lt_" << i << ' ' << vertices[i].q << ' ' << vertices[i].Q << ' '
             << vertices[i].coords[0] << ' ' << vertices[i].coords[1] << ' '
             << vertices[i].coords[2] << ' ' << vertices[i].Phi << '\n';
    }
    fout.close();
    fout.open(path_to_data / "edge_table.txt");
    fout << "id" << ' ' << "from" << ' ' << "to" << ' ' << "current" << ' ' << "sigma" << '\n';
    for (int i = 0; i < edges.size(); i++) {
        if (!edges_activity[i])
            continue;
        fout << "lt_" << i << ' ' << "lt_" << edges[i].from << ' ' << "lt_" << edges[i].to << ' '
             << edges[i].current << ' ' << edges[i].sigma << '\n';
    }
    fout.close();
}

void LightningTree::ReturnPhi(const std::filesystem::path& path_to_data,
                              const std::array<double, 3>& start_r,
                              const std::array<double, 3>& end_r) {
    std::ofstream fout(path_to_data / "phi_info.txt");
    fout << 'z' << ' ' << "full_phi" << ' ' << "ext_phi" << '\n';
    for (double z = start_r[2]; z < end_r[2]; z += h) {
        double full_sum = 0;
        double ext_sum = 0;
        for (double x = start_r[0]; x < end_r[0]; x += h) {
            for (double y = start_r[1]; y < end_r[1]; y += h) {
                full_sum += Potential({x, y, z}) +
                            external_field_potential({x, y, z});  // Стоит ли считать потенциал в
                                                                  // каждой точке плоскости
                                                                  // относительно всех зарядов???
                ext_sum += external_field_potential({x, y, z});
            }
        }
        fout << z << ' ' << full_sum * (h * h) / ((end_r[0] - start_r[0]) * (end_r[1] - start_r[1]))
             << ' ' << ext_sum * (h * h) / ((end_r[0] - start_r[0]) * (end_r[1] - start_r[1]))
             << '\n';
    }
    fout.close();
}
