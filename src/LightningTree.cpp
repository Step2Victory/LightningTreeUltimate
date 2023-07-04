#include "LightningTree.h"

LightningTree::LightningTree(
    double h, double delta_t, double r, double R,
    size_t periphery_size, double q_plus_max, double q_minus_max,
    double Q_plus_s, double Q_minus_s, double resistance,
    double E_plus, double E_minus, double alpha, double beta,
    double sigma,
    std::function<double(double, double, double)>
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
    for(auto& point : vertices){
        double Phi = 0;
        for (auto& vertex : vertices)
        {
            double l = abs(sqrt(pow((vertex.coords[0] - point.coords[0]), 2) + pow((vertex.coords[1] - point.coords[1]), 2) + pow((vertex.coords[2] - point.coords[2]), 2)));
            double mirror_l = abs(sqrt(pow((vertex.coords[0] - point.coords[0]), 2) + pow((vertex.coords[1] - point.coords[1]), 2) + pow((vertex.coords[2] + point.coords[2]), 2)));
            if (l < kEps)
            {
                Phi += vertex.q / (4 * PI * epsilon_0) * (1 / (h + r) + 1 / (mirror_l + r)) + 
                            vertex.Q / (4 * PI * epsilon_0) * (1 / (h + R) + 1 / (mirror_l + R));
            }
            Phi += vertex.q / (4 * PI * epsilon_0) * (1 / (l + r) + 1 / (mirror_l + r)) +
                        vertex.Q / (4 * PI * epsilon_0) * (1 / (l + R) + 1 / (mirror_l + R));
        }
        // ???
        point.Phi = Phi + phi_a->getValue(point.coords); 
    }
}

double LightningTree::CountElectricity(const size_t v_from_id, const size_t v_to_id) const {
    /*
    Описание метода
    */
    double l = abs(sqrt(pow((vertices[v_from_id].coords[0] - vertices[v_to_id].coords[0]), 2) + 
                        pow((vertices[v_from_id].coords[1] - vertices[v_to_id].coords[1]), 2) + 
                        pow((vertices[v_from_id].coords[2] - vertices[v_to_id].coords[2]), 2)));
    double phi_from = vertices[v_from_id].Phi;
    double phi_to = vertices[v_to_id].Phi;
    // ???
    if (std::isinf((phi_from - phi_to) / l))
    {
        throw std::runtime_error{ "Electric field along edge is infinity!" };
    }
    return (phi_from - phi_to) / l;
}

void LightningTree::CountSigma() {
    /*
    Описание метода
    */
    for (auto edge : edges) {
        double E = CountElectricity(edge.from, edge.to);
        edge.sigma = edge.sigma * std::exp((alpha * E * E - beta) * delta_t);
        // ???
        if (std::isinf(sigma))
        {
            throw std::runtime_error{ "Sigma is infinity!" };
        }
    }
}

void LightningTree::CountCurrent() {
    /*
    Описание метода
    */
    for (auto edge : edges) {
        edge.current = PI * r * r * edge.sigma * CountElectricity(edge.from, edge.to);
    }
}

void LightningTree::countCoords(std::array<double, 3>& result, const size_t vertex_id, const std::vector<int> point) {
    result = {
        vertices[vertex_id].coords[0] + (1 - point[0]) * h,
        vertices[vertex_id].coords[1] + (1 - point[1]) * h,
        vertices[vertex_id].coords[2] + (1 - point[2]) * h
    };
}

// ???
cubic_grid LightningTree::CreateNode(size_t vertex, size_t edge, const std::vector<int>& point) {
    /*
    Описание метода
    */
   // TO DO 
   cubic_grid node;
   for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                node[i][j][k] = -1;
                
            }
        }
    }
    node[1][1][1] = vertex;
    node[2-point[0]][2-point[1]][2-point[2]] = edge;
    return node;
}

void LightningTree::Transport() {
    /*
    Описание метода. Узнать что делает этот метод!
    */
    // TO DO
    std::unordered_map<size_t, std::pair<double, double>> delta_charges;
    for (auto& elem : graph)
    {
        size_t vertex_id = elem[1][1][1];
        // std::vector<EdgePtr> edges = elem.second;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    if (elem[i][j][k] == -1 && (i == 1 && j == 1 && k == 1)) continue;
                    if (edges[elem[i][j][k]].from == vertex_id)
                    {
                        delta_charges[vertex_id].first -= edges[elem[i][j][k]].current;
                    }
                    else
                    {
                        delta_charges[vertex_id].first += edges[elem[i][j][k]].current;
                    }
                }
            }
        }
        // delta_charges[vertex_id].first -= CurrentSheath(vertex_id);
        // delta_charges[vertex_id].second += CurrentSheath(vertex_id);
    }
    // for (auto& elem : delta_charges)
    // {
    //     vertices[elem.first].q += elem.second.first * delta_t;
    //     vertices[elem.first].Q += elem.second.second * delta_t;
    // }
    iter_number++;
}

void LightningTree::Grow() {
    /*
    Описание метода
    */
    // TO DO
    std::vector<cubic_grid> temp_graph;
    for (unsigned v = graph.size(); v-- > 0; ) {
        size_t vertex_id = graph[v][1][1][1];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {

                    if (graph[v][i][j][k] != -1 && (i == 1 && j == 1 && k == 1)) continue;

                    std::array<double, 3> coords = {0, 0, 0};
                    countCoords(coords, vertex_id, std::vector{ i, j, k });
                    LightningTree::Vertex new_vertex = { 0, 0, 0, coords, 0 };
                    vertices.push_back(new_vertex);

                    if (GrowthCriterion(vertex_id, vertices.size() - 1)) {
                        LightningTree::Edge new_edge = { vertex_id, vertices.size() - 1, 0, 0};
                        edges.push_back(new_edge);
                        graph[v][i][j][k] = edges.size() - 1;
                        temp_graph.push_back(CreateNode(vertices.size() - 1, edges.size() - 1, std::vector{i, j, k}));
                    }
                    else {
                        vertices.pop_back();
                    }
                }
            }
        }
    }
    graph.insert(graph.end(), temp_graph.begin(), temp_graph.end());
}

void LightningTree::Delete() {
    /*
    Описание метода
    */
    // TO DO

}

bool LightningTree::GrowthCriterion(const size_t v_from_id, const size_t v_to_id) const {
    // ???
    double probability = dis(gen);
    double E = CountElectricity(v_from_id, v_to_id);
    if (E > E_plus)
    {
        return (1 - std::exp(-std::pow(((E - E_plus) / E_plus), 2.5))) > probability;
    }
    else if (-E > E_minus)
    {
        return (1 - std::exp(-std::pow(((-E - E_minus) / E_minus), 2.5))) > probability;
    }
    return false;
}
bool LightningTree::DeletionCriterion(size_t vertex_id) const {
    /*
    Описание метода
    */
    // TO DO
}