#include "LightningTree.h"

LightningTree::LightningTree(
    double _h, double _delta_t, double _r, double _R,
    size_t _periphery_size, double _q_plus_max, double _q_minus_max,
    double _Q_plus_s, double _Q_minus_s, double _resistance,
    double _E_plus, double _E_minus, double _alpha, double _beta,
    double _sigma,
    std::function<double(double, double, double)>
        _external_field_potential)
    : h(_h),
      delta_t(_delta_t),
      r(_r),
      R(_R),
      periphery_size(_periphery_size),
      q_plus_max(_q_plus_max),
      q_minus_max(_q_minus_max),
      Q_plus_s(_Q_plus_s),
      Q_minus_s(_Q_minus_s),
      resistance(_resistance),
      E_plus(_E_plus),
      E_minus(_E_minus),
      alpha(_alpha),
      beta(_beta),
      sigma(_sigma),
      external_field_potential(_external_field_potential) {
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
                Phi = vertex.q / (4 * M_PI * epsilon_0) * (1 / (h + r) + 1 / (h + mirror_l + r)) + 
                            vertex.Q / (4 * M_PI * epsilon_0) * (1 / (h + R) + 1 / (h + mirror_l + R));
            }
            Phi = vertex.q / (4 * M_PI * epsilon_0) * (1 / (l + r) + 1 / (mirror_l + r)) +
                        vertex.Q / (4 * M_PI * epsilon_0) * (1 / (l + R) + 1 / (mirror_l + R));
        }
        point.Phi = Phi + phi_a->getValue(point.coords);
    }
}

double LightningTree::CountElectricity(const Edge& edge) {
    /*
    Описание метода
    */
    double l = abs(sqrt(pow((vertisec[edge.from].coords[0] - vertisec[edge.to].coords[0]), 2) + 
                        pow((vertisec[edge.from].coords[1] - vertisec[edge.to].coords[1]), 2) + 
                        pow((vertisec[edge.from].coords[2] - vertisec[edge.to].coords[2]), 2)));
    double phi_from = vertisec[edge.from].Phi;
    double phi_to = vertisec[edge.to].Phi;
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
        double E = CountElectricity(edge);
        edge.sigma = edge.sigma * std::exp((alpha * E * E - beta) * delta_t);
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
        edge.current = M_PI * r * r * edge.sigma * CountElectricity(edge);
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

void LightningTree::countCoords(std::array<double, 3>& result, size_t vertex_id, const std::vector<int> point) {
    result = {
        vertices[vertex_id].coords[0] + (1 - point[0]) * h,
        vertices[vertex_id].coords[1] + (1 - point[1]) * h,
        vertices[vertex_id].coords[2] + (1 - point[2]) * h
    }
}

void LightningTree::Grow() {
    /*
    Описание метода
    */
    // TO DO
    for (unsigned v = graph.size(); v-- > 0; ) {
        size_t vertex_id = graph[v][1][1][1];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    if (graph[v][i][j][k] != -1 && (i == 1 && j == 1 && k == 1)) continue;
                    std::array<double, 3> coords = {0, 0, 0}
                    countCoords(coords, vertex_id, std::vector{ i, j, k });
                    LightningTree::Vertex new_vertex = { 0, 0, 0, coords, 0 };
                    vertices.push_back(new_vertex);
                    LightningTree::Edge new_edge = { vertex_id, vertices.size() - 1, 0, 0, 0 };
                    if (GrowthCriterion(new_edge)) {
                        edges.push_back(new_edge);
                        graph[v][i][j][k] = edges.size() - 1;
                        graph.push_back(CreateNode(vertices.size() - 1), edges.size() - 1, std::vector{i, j, k});
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
    double E = CountElectricity(edge);
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

LightningTree::cubic_grid LightningTree::CreateNode(size_t vertex, size_t edge, const std::vector<int>& point) {
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