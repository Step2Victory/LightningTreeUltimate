#include "LightningTree.h"

#include "yaml-cpp/yaml.h"
#include <glog/logging.h>

#include <numbers>
#include <cmath>
#include <fstream>
#include <string>
#include <exception>
#include <optional>

#include "ExternalField.h"
#include "RandomDirections.h"

bool checkDouble(double value) {
    if (std::isnan(value) || std::isinf(value)) {
        return false;
    }
    return true;
}

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
    // std::cout<<"Расчёт внешнего поля\n";
    // auto start = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    // external_field_potential = countExternalField(layers, start_r, end_r, h);
    external_field_potential = countExternalField_octree(layers, start_r, end_r, h); // версия с октодеревом
    // external_field_potential = constExternalField(6000.0);
    // auto end = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    // std::cout<<"Расчёт завершён. Время = "<<end-start<<'\n';
    // std::cout << "external field max: " << external_field_potential({0, 0, 10000}) << '\n';

    seed = config["seed"].as<int>();
    degree_probability_growth = config["degree_probability_growth"].as<double>();
    gen = std::mt19937(seed);
    dis = std::uniform_real_distribution<>(0, 1);
    iter_number = 0;
    max_number_edges = config["max_number_edges"].as<int>();
    // dynoctree = DynamicOctree(start_r, end_r);
    auto first =
        addVertex(Vertex{.q = 0,
                         .Q = 0,
                         .Phi = 0,
                         .coords = {(end_r[0] + start_r[0]) / 2 + h/2, (end_r[1] + start_r[1]) / 2,
                                    (end_r[2] + start_r[2]) / 2 + h/2},
                         .internal_coords = {0, 0, 0},
                         .number_edges = 0,
                         .growless_iter_number = 0});
    // dynoctree.add_charge(first, vertices);
    auto second =
        addVertex(Vertex{.q = 0,
                         .Q = 0,
                         .Phi = 0,
                         .coords = {(end_r[0] + start_r[0]) / 2 - h/2, (end_r[1] + start_r[1]) / 2,
                                    (end_r[2] + start_r[2]) / 2 - h/2},
                         .internal_coords = {0, 0, 1},
                         .number_edges = 0,
                         .growless_iter_number = 0});
    // dynoctree.add_charge(second, vertices);
    addEdge(first, second);

    // // dynoctree.print('d');
}

LightningTree::LightningTree(double h_, double delta_t_, double r_, double R_, size_t periphery_size_,
                             double q_plus_max_, double q_minus_max_, double Q_plus_s_,
                             double Q_minus_s_, double resistance_, double E_plus_, double E_minus_,
                             double alpha_, double beta_, double sigma_, std::array<double, 3> start_r_,
                             std::array<double, 3> end_r_, double degree_probability_growth_,
                             int seed_, int max_number_edges_)
    : start_r(start_r_),
      end_r(end_r_),
      h(h_),
      delta_t(delta_t_),
      r(r_),
      R(R_),
      q_plus_max(q_plus_max_),
      q_minus_max(q_minus_max_),
      Q_plus_s(Q_plus_s_),
      Q_minus_s(Q_minus_s_),
      resistance(resistance_),
      E_plus(E_plus_),
      E_minus(E_minus_),
      alpha(alpha_),
      beta(beta_),
      sigma(sigma_),
      //   external_field_potential(constExternalField(100000.0)),
      degree_probability_growth(degree_probability_growth_),
      periphery_size(periphery_size_),
      external_field_potential(countExternalField(
          std::vector({ChargeLayer{.p_0 = -0.0000007, .h = h, .L = 200, .r = start_r, .a = 1},
                       ChargeLayer{.p_0 = 0.0000007, .h = h, .L = 200, .r = end_r, .a = 1}}),
          std::array<double, 3>{end_r[0] - start_r[0] - 5 * h, end_r[1] - start_r[1] - 5 * h,
                                end_r[2] - start_r[2] - 20 * h},
          std::array<double, 3>{end_r[0] - start_r[0] + 5 * h, end_r[1] - start_r[1] + 5 * h,
                                end_r[2] - start_r[2] + 20 * h},
          h)),
      seed(seed_),
      max_number_edges(max_number_edges_) {
    gen = std::mt19937(seed_);
    dis = std::uniform_real_distribution<>(0, 1);
    iter_number = 0;
    // dynoctree = DynamicOctree(start_r, end_r);
    auto first =
        addVertex(Vertex{.q = 0,
                         .Q = 0,
                         .Phi = 0,
                         .coords = {(end_r[0] + start_r[0]) / 2, (end_r[1] + start_r[1]) / 2,
                                    (end_r[2] + start_r[2]) / 2 - h},
                         .internal_coords = {0, 0, 0},
                         .number_edges = 0,
                         .growless_iter_number = 0});
    // dynoctree.add_charge(first, vertices);
    auto second =
        addVertex(Vertex{.q = 0,
                         .Q = 0,
                         .Phi = 0,
                         .coords = {(end_r[0] + start_r[0]) / 2, (end_r[1] + start_r[1]) / 2,
                                    (end_r[2] + start_r[2]) / 2},
                         .internal_coords = {0, 0, 1},
                         .number_edges = 0,
                         .growless_iter_number = 0});
    // dynoctree.add_charge(second, vertices);
    addEdge(first, second);
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
        double k = 1 / (4 * std::numbers::pi * epsilon_0);
        if (l < kEps) {
            if (!checkDouble(vertex.q) && checkDouble(vertex.Q)) {
                LOG(INFO) << "Q = " << vertex.Q << ", q = " << vertex.q;
            }
            Phi += vertex.q * k * (1 / (h / 2 + r) - 1 / (mirror_l + r)) +
                   vertex.Q * k * (1 / (h / 2 + R) - 1 / (mirror_l + R));
        } else {
            if (!checkDouble(vertex.q) && checkDouble(vertex.Q)) {
                LOG(INFO) << "Q = " << vertex.Q << ", q = " << vertex.q;
            }
            Phi += vertex.q * k * (1 / (l + r) - 1 / (mirror_l + r)) +
                   vertex.Q * k * (1 / (l + R) - 1 / (mirror_l + R));
        }
    }
    // Phi = dynoctree.potencial_in_point(coords, r, R, h);
    // std::cout<<"Потенциал = "<<Phi<<"\n";
    return Phi;
}

void LightningTree::CountPotential() {
    for (auto& point : vertices) {
        point.Phi = Potential(point.coords) + external_field_potential(point.coords);
    }
}

double LightningTree::CountElectricity(const Vertex& vertex_from, const Vertex& vertex_to) const {
    double l = countDistance(vertex_from.coords, vertex_to.coords);
    double phi_from = vertex_from.Phi;
    double phi_to = vertex_to.Phi;

    if (!checkDouble(phi_from - phi_to) / l) {
        LOG(INFO) << "Electric field in from "
                  << ": " << phi_from << ". "
                  << "Electric field in to"
                  << ": " << phi_to;
        throw std::runtime_error{"Electric field along edge is inf or nan!"};
    }

    return -(phi_from - phi_to) / l;
}

void LightningTree::CountSigma() {
    for (auto& edge : edges) {
        double E = CountElectricity(vertices[edge.from], vertices[edge.to]);
        // std::cout << E << "\n";
        edge.sigma *= std::exp((alpha * E * E - beta) * delta_t);
        // std::cout<<"alpha = "<<alpha<<" beta = "<<beta<<" E = "<<E<<" sima = "<<edge.sigma<<std::endl;

        if (!checkDouble(sigma)) {
            LOG(INFO) << "Incorrect sigma value " << sigma;
            throw std::runtime_error{"Sigma in edge is inf or nan!"};
        }
    }
}

void LightningTree::CountCurrent() {
    for (auto& edge : edges) {
        double E = CountElectricity(vertices[edge.from], vertices[edge.to]);
        edge.current = std::numbers::pi * r * r * edge.sigma * E;
                    //    CountElectricity(vertices[edge.from], vertices[edge.to]);
        // std::cout<<"Радиус канала = "<<r<<" Сигма = "<<edge.sigma<<" Напряженность = "<<E<<" Ток = "<<edge.current<<std::endl;
        if (!checkDouble(edge.current)) {
            LOG(INFO) << "Error in current " << edge.current;
            throw std::runtime_error{"bad current"};
        }
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

// double LightningTree::countDistance(const std::array<double, 3>& r_point,
//                                     const std::array<double, 3>& l_point) const {
//     return std::abs(std::sqrt(std::pow((r_point[0] - l_point[0]), 2) +
//                               std::pow((r_point[1] - l_point[1]), 2) +
//                               std::pow((r_point[2] - l_point[2]), 2)));
// }

LightningTree::cubic_grid LightningTree::CreateEmptyNode() {
    cubic_grid node;
    node[1][1][1] = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                node[i][j][k] = -1;
            }
        }
    }
    node[1][1][1] = 0;
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
    // for(size_t i = 0; i < vertices.size(); i++){
    //     if(vertices_activity[i])
    //         std::cout<<"q = "<<vertices[i].q<<", Q = "<<vertices[i].Q<<'\n';
    //     else 
    //         std::cout<<"Вершина "<<i<<" не активна\n";
    // }
    iter_number++;
    // dynoctree.recalc_sumCharge(vertices); // Актуализация зарядов в октодереве
    // // dynoctree.print('d');
}

void LightningTree::Grow() {
    std::vector<cubic_grid> temp_graph;
    int number_of_vertices_before_grow = graph.size();
    for (size_t v_from_id = number_of_vertices_before_grow; v_from_id-- > 0;) {
        // std::cout << v_from_id << "\t" << vertices_activity.size() << std::endl;
        if (!vertices_activity[v_from_id] || vertices[v_from_id].number_edges >= static_cast<size_t>(max_number_edges))
            continue;
        bool notGrow = true;
        // auto directions = randomDirections();
        auto temp = GetPotencialDir(v_from_id);
        auto directions = sortDir(temp);
        // for(size_t i = 0; i < temp.size(); i++){
        //     std::cout<<directions[i].first<<", ";
        // }
        // std::cout<<std::endl;
        for (const auto& dir : directions) {
            if (vertices[v_from_id].number_edges >= static_cast<size_t>(max_number_edges))
                break;
            // int i = dir[0];
            // int j = dir[1];
            // int k = dir[2];
            int i = dir.second[0];
            int j = dir.second[1];
            int k = dir.second[2];
            int e_id = graph[v_from_id][i][j][k];

            if ((i == 1 && j == 1 && k == 1) || (e_id != -1 && edges_activity[e_id]))
                continue;

            auto internal_coords = countInternalCoords(v_from_id, dir.second);
            auto coords = countCoords(v_from_id, dir.second);
            std::optional<size_t> v_to_id;
            Vertex vertex;

            if (internal_coords_to_id.contains(internal_coords)) {
                v_to_id = internal_coords_to_id[internal_coords];
                vertex = vertices[*v_to_id];
            } else {
                vertex = Vertex{.q = 0,
                                .Q = 0,
                                .Phi = dir.first,
                                // .Phi = Potential(coords),
                                .coords = coords,
                                .internal_coords = internal_coords,
                                .number_edges = 0,
                                .growless_iter_number = 0};
            }
            if (v_to_id.has_value() && vertices_activity[*v_to_id]) {
                continue;
            }
            bool criterion_value = GrowthCriterion(vertices[v_from_id], vertex);
            if (criterion_value) {
                vertices[v_from_id].number_edges++;
                if (!v_to_id.has_value()) {
                    v_to_id = addVertex(vertex);
                } else {
                    vertices[*v_to_id].number_edges++;
                    vertices_activity[*v_to_id] = true;
                }

                // dynoctree.add_charge(*v_to_id, vertices);

                if (e_id == -1) {
                    addEdge(v_from_id, *v_to_id);
                } else {
                    edges_activity[e_id] = true;
                }
                notGrow = false;
            }
        }
        if (notGrow) {
            vertices[v_from_id].growless_iter_number++;
        }
    }
    graph.insert(graph.end(), temp_graph.begin(), temp_graph.end());
}

size_t LightningTree::addVertex(Vertex vertex) {
    vertex.number_edges++;
    vertices.push_back(vertex);
    internal_coords_to_id[vertex.internal_coords] = vertices.size() - 1;
    vertices_activity.push_back(true);
    vertices_peripherality.push_back(true);
    graph.push_back(CreateEmptyNode());
    return vertices.size() - 1;
}

size_t LightningTree::addEdge(size_t from, size_t to) {
    edges.push_back({.from = from, .to = to, .current = 0, .sigma = sigma});
    edges_activity.push_back(true);
    std::array<int, 3> dir = {
        1 + (vertices[to].internal_coords[0] - vertices[from].internal_coords[0]),
        1 + (vertices[to].internal_coords[1] - vertices[from].internal_coords[1]),
        1 + (vertices[to].internal_coords[2] - vertices[from].internal_coords[2])};
    graph[from][dir[0]][dir[1]][dir[2]] = edges.size() - 1;
    graph[to][2 - dir[0]][2 - dir[1]][2 - dir[2]] = edges.size() - 1;
    return edges.size() - 1;
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
            // dynoctree.delete_charge(v_id, vertices);
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
                            vertices_activity[edges[e_id].from] = false;
                        }
                        // if (vertices[edges[e_id].to].number_edges == 1) {
                        //     vertices_peripherality[edges[e_id].to] = true;
                        // }
                        // if (vertices[edges[e_id].to].number_edges == 0) {
                        //     vertices_activity[edges[e_id].to] = false;
                        // }
                    }
                }
            }
        }
    }
}

bool LightningTree::GrowthCriterion(const Vertex& vertex_from, const Vertex& vertex_to) const {
    /*
    Описание метода
    */
    double probability = dis(gen);
    double E = CountElectricity(vertex_from, vertex_to);
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
    for (size_t i = 0; i < vertices.size(); i++) {
        fout << "lt_" << i << ' ' << vertices[i].q << ' ' << vertices[i].Q << ' '
             << vertices[i].coords[0] << ' ' << vertices[i].coords[1] << ' '
             << vertices[i].coords[2] << ' ' << vertices[i].Phi << '\n';
    }
    fout.close();
    fout.open(path_to_data / "edge_table.txt");
    fout << "id" << ' ' << "from" << ' ' << "to" << ' ' << "current" << ' ' << "sigma" << '\n';
    for (size_t i = 0; i < edges.size(); i++) {
        if (!edges_activity[i])
            continue;
        fout << "lt_" << i << ' ' << "lt_" << edges[i].from << ' ' << "lt_" << edges[i].to << ' '
             << edges[i].current << ' ' << edges[i].sigma << '\n';
    }
    fout.close();
}

void LightningTree::ReturnPhi(const std::filesystem::path& path_to_data) {
    std::ofstream fout(path_to_data / "phi_info.txt");
    fout << 'z' << ' ' << "full_phi" << ' ' << "ext_phi" << '\n';
    for (double z = start_r[2]; z < end_r[2]; z += h) {
        double full_sum = 0;
        double ext_sum = 0;
        for (double x = start_r[0]; x < end_r[0]; x += h) {
            for (double y = start_r[1]; y < end_r[1]; y += h) {
                full_sum += Potential({x, y, z}) + external_field_potential({x, y, z});
                ext_sum += external_field_potential({x, y, z});
            }
        }
        fout << z << ' ' << full_sum * (h * h) / ((end_r[0] - start_r[0]) * (end_r[1] - start_r[1]))
             << ' ' << ext_sum * (h * h) / ((end_r[0] - start_r[0]) * (end_r[1] - start_r[1]))
             << '\n';
    }
    fout.close();
}

std::array<std::pair<double, std::array<int, 3>>, 26> LightningTree::GetPotencialDir(size_t vid)
{
    std::array<std::pair<double, std::array<int, 3>>, 26> potencial_dir;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            for (int k = 0; k < 3; ++k)
            {
                if (i == 1 && j == 1 && k == 1)
                {
                    continue;
                }
                std::array<double, 3> new_point = countCoords(vid, {i, j, k});
                double l = countDistance(vertices[vid].coords, new_point);
                double E = (-(vertices[vid].Phi - Potential(new_point)) / l);
                // std::cout<<Phi<<", ";
                potencial_dir[i*9 + j*3 + k] = std::pair<double, std::array<int, 3>>(E, {i, j, k});
            }
        }
    }
    // std::cout<<'\n';
    return potencial_dir;
}