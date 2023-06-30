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
    CountElectricity();
    CountSigma();
    CountCurrent();
    double grow_time_period = 0.001;
    int grow_iter_period = grow_time_period / delta_t;
    if (iter_number % grow_iter_period == 0) {
        Grow();
        Delete();
    }
}

void LightningTree::CountSigma() {
    /*
    Описание метода
    */
    for (auto edge : edges) {
        CountElectricity(edge);
        edge.sigma = edge.sigma * std::exp((alpha * edge.E * edge.E - beta) * delta_t);
        if (std::isinf(sigma))
        {
            throw std::runtime_error{ "Sigma is infinity!" };
        }
    }
}

double LightningTree::CountPotential_q(std::array<double, 3> point) {
    /*
    * Описание метода
    */
    for (auto& vertex : vertices)
    {
        double l = abs(sqrt(pow((vertex.coords[0] - point[0]), 2) + pow((vertex.coords[1] - point[1]), 2) + pow((vertex.coords[2] - point[2]), 2)));
        double mirror_l = abs(sqrt(pow((vertex.coords[0] - point[0]), 2) + pow((vertex.coords[1] - point[1]), 2) + pow((vertex.coords[2] + point[2]), 2)));
        if (l < kEps)
        {
            return vertex.q / (4 * M_PI * epsilon_0) * (1 / (h + r) + 1 / (h + mirror_l + r));
        }
        return vertex.q / (4 * M_PI * epsilon_0) * (1 / (l + r) + 1 / (mirror_l + r));
    }
}

double LightningTree::CountPotential_Q(std::array<double, 3> point) {
    /*
    * Описание метода
    */
    for (auto& vertex : vertices)
    {
        double L = abs(sqrt(pow((vertex.coords[0] - point[0]), 2) + pow((vertex.coords[1] - point[1]), 2) + pow((vertex.coords[2] - point[2]), 2)));
        double mirror_L = abs(sqrt(pow((vertex.coords[0] - point[0]), 2) + pow((vertex.coords[1] - point[1]), 2) + pow((vertex.coords[2] + point[2]), 2)));
        if (L < kEps)
        {
            return vertex.Q / (4 * M_PI * epsilon_0) * (1 / (h + R) + 1 / (h + mirror_L + R));
        }
        return vertex.Q / (4 * M_PI * epsilon_0) * (1 / (L + R) + 1 / (mirror_L + R));
    }
}

void LightningTree::CountElectricity() {
    /*
    Описание метода
    */
    for (auto edge : edges) {
        double l = abs(sqrt(pow((vertices[edge.from].coords[0] - vertices[edge.to].coords[0]), 2) + pow((vertices[edge.from].coords[1] - vertices[edge.to].coords[1]), 2) + pow((vertices[edge.from].coords[2] - vertices[edge.to].coords[2]), 2)));
        // double l = h;
        double phi_from = CountPotential_q(vertices[edge.from].coords) + CountPotential_Q(vertices[edge.from].coords) + phi_a->getValue(vertices[edge.from].coords); // формула (6)
        double phi_to = CountPotential_q(vertices[edge.to].coords) + CountPotential_Q(vertices[edge.to].coords) + phi_a->getValue(vertices[edge.to].coords);
        if (std::isinf((phi_from - phi_to) / l))
        {
            throw std::runtime_error{ "Electric field along edge is infinity!" };
        }
        edge.E = (phi_from - phi_to) / l;
    }
}

void LightningTree::CountElectricity(Edge edge) {
    /*
    Описание метода
    */
    double l = abs(sqrt(pow((vertices[edge.from].coords[0] - vertices[edge.to].coords[0]), 2) + pow((vertices[edge.from].coords[1] - vertices[edge.to].coords[1]), 2) + pow((vertices[edge.from].coords[2] - vertices[edge.to].coords[2]), 2)));
    // double l = h;
    double phi_from = CountPotential_q(vertices[edge.from].coords) + CountPotential_Q(vertices[edge.from].coords) + phi_a->getValue(vertices[edge.from].coords); // формула (6)
    double phi_to = CountPotential_q(vertices[edge.to].coords) + CountPotential_Q(vertices[edge.to].coords) + phi_a->getValue(vertices[edge.to].coords);
    if (std::isinf((phi_from - phi_to) / l))
    {
        throw std::runtime_error{ "Electric field along edge is infinity!" };
    }
    edge.E = (phi_from - phi_to)/l;
}

//double LightningTree::CountElectricity(size_t vertex_id, const std::array<double, 3> coords) {
//    /*
//    Описание метода
//    */
//    double l = abs(sqrt(pow((vertices[vertex_id].coords[0] - coords[0]), 2) + 
//                        pow((vertices[vertex_id].coords[1] - coords[1]), 2) + 
//                        pow((vertices[vertex_id].coords[2] - coords[2]), 2)));
//    // double l = h;
//    double phi_from = CountPotential_q(vertices[vertex_id].coords) + CountPotential_Q(vertices[vertex_id].coords) + phi_a->getValue(vertices[vertex_id].coords); // формула (6)
//    double phi_to = CountPotential_q(coords) + CountPotential_Q(coords) + phi_a->getValue(coords);
//    if (std::isinf((phi_from - phi_to) / l))
//    {
//        throw std::runtime_error{ "Electric field along edge is infinity!" };
//    }
//    return (phi_from - phi_to) / l;
//}

void LightningTree::CountCurrent() {
    /*
    Описание метода
    */
    for (auto edge : edges) {
        edge.current = M_PI * r * r * edge.sigma * edge.E;
    }
}

void LightningTree::Transport() {
    /*
    Описание метода. Узнать что делает этот метод!
    */
    // TO DO
    std::unordered_map<VertexPtr, std::pair<double, double>> delta_charges;
    for (auto& elem : graph)
    {
        VertexPtr vertex = elem.first;
        std::vector<EdgePtr> edges = elem.second;
        for (auto& edge : edges)
        {
            if (edge->from == vertex)
            {
                delta_charges[vertex].first -= edge.current;
            }
            else
            {
                delta_charges[vertex].first += edge.current;
            }
        }
        delta_charges[vertex].first -= CurrentSheath(vertex);
        delta_charges[vertex].second += CurrentSheath(vertex);
    }
    for (auto& elem : delta_charges)
    {
        elem.first->q += delta_charges[elem.first].first * delta_t;
        elem.first->Q += delta_charges[elem.first].second * delta_t;
    }
    iter_number++;
}

std::array<double, 3> countCoords(size_t vertex_id, std::vector<int> point) {
    std::array<double, 3> coords = {
        vertices[vertex_id].coords[0] + (1 - point[0]) * h,
        vertices[vertex_id].coords[1] + (1 - point[1]) * h,
        vertices[vertex_id].coords[2] + (1 - point[2]) * h
    }
    return coords;
}

void LightningTree::Grow() {
    /*
    Описание метода
    */
    // TO DO
    for (unsigned v = graph.size(); v-- > 0; ) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    if (graph[v][i][j][k] != -1 && (i == 1 && j == 1 && k == 1)) continue;
                    std::array<double, 3> coords = countCoords(graph[v][1][1][1], std::vector{ i, j, k });
                    Vertex new_vertex = { 0, 0, coords, 0 };
                    vertices.push_back(new_vertex);
                    Edge new_edge = { graph[v][i][j][k], vertices.size() - 1, 0, 0, 0 };
                    if (GrowthCriterion(new_edge)) {
                        edges.push_back(new_edge);
                        graph[v][i][j][k] = vertices.size() - 1;
                        graph.push_back(CreateNode(vertices.size() - 1));
                    }
                    else {
                        vertices.pop_back();
                    }
                }
            }
        }
    }
}

void LightningTree::Delete() {
    /*
    Описание метода
    */
    // TO DO
}

bool LightningTree::GrowthCriterion(Edge edge) const {
    double probability = dis(gen);
    edge.E = CountElectricity(edge);
    if (edge.E > E_plus)
    {
        return (1 - std::exp(-std::pow(((edge.E - E_plus) / E_plus), 1))) > probability;
    }
    else if (-edge.E > E_minus)
    {
        return (1 - std::exp(-std::pow(((-edge.E - E_minus) / E_minus), 1))) > probability;
    }
    return false;
}

bool LightningTree::GrowthCriterion(size_t vertex_id, const std::array<double, 3>& coords) const {
    /*
    Описание метода
    */
    double probability = dis(gen);
    double E = CountElectricity(vertex_id, coords);
    if (E > E_plus)
    {
        return (1 - std::exp(-std::pow(((E - E_plus) / E_plus), 1))) > probability;
    }
    else if (-E > E_minus)
    {
        return (1 - std::exp(-std::pow(((-E - E_minus) / E_minus), 1))) > probability;
    }
    return false;
}

bool LightningTree::DeletionCriterion(size_t vertex_id) const {
    /*
    Описание метода
    */
    // TO DO
}

LightningTree::cubic_grid LightningTree::CreateNode(size_t vertex) {
    // TO DO 
}