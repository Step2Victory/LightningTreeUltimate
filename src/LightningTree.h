#pragma once
#include "DynamicOctree.h"

#include <vector>
#include <array>
#include <functional>
#include <filesystem>
#include <random>
#include <unordered_map>
#include <map>

class LightningTree {
    struct Vertex;

private:
    using cubic_grid = std::array<std::array<std::array<int, 3>, 3>, 3>;

public:
    LightningTree(const std::filesystem::path& path_to_config_file);
    LightningTree(double h, double delta_t, double r, double R, size_t periphery_size,
                  double q_plus_max, double q_minus_max, double Q_plus_s, double Q_minus_s,
                  double resistance, double E_plus, double E_minus, double alpha, double beta,
                  double sigma, std::array<double, 3> start_r, std::array<double, 3> end_r,
                  double degree_probability_growth, int seed, int max_number_edges);
    void NextIter();
    void CountSigma();
    double Potential(const std::array<double, 3>&);
    void CountPotential();
    double CountElectricity(const Vertex& from, const Vertex& to) const;
    void CountCurrent();
    std::array<double, 3> countCoords(size_t, const std::array<int, 3>&);
    // std::array<int, 3> countInternalCoords(size_t, const std::array<int, 3>&);
    std::array<int, 3> count_internal_coords(const std::array<double, 3>&);
    size_t get_hash_id(const std::array<int, 3>&);

    // double countDistance(const std::array<double, 3>&, const std::array<double, 3>&) const;
    cubic_grid CreateNode(size_t, const std::array<int, 3>&);
    cubic_grid CreateEmptyNode();
    size_t addVertex(Vertex);
    size_t addEdge(size_t, size_t);

    void Transport();
    void Grow();
    void Delete();

    bool GrowthCriterion(const Vertex& vertex_from, const Vertex& vertex_to) const;
    bool DeletionCriterion(size_t) const;

    void AllParams() const;
    void Info() const;
    void WriteResponse(int response) const;
    void ReturnFiles(const std::filesystem::path&);
    void ReturnPhi(const std::filesystem::path&);

    std::array<std::pair<double, std::array<int, 3>>, 26> GetPotencialDir(size_t);


    
private:
    bool TryAddEdge(size_t v_from_id, const std::array<int, 3>& dir);

    struct Vertex {
        // структура вершины
        double q; // заряд в точке
        double Q; // заряд в чехле
        double Phi; // потенциал в вершине
        std::array<double, 3> coords; // координаты вершины
        std::array<int, 3> internal_coords; // относительные координаты вершины
        size_t number_edges; // количетсво рёбер
        int growless_iter_number; // Количество пропущенных итераций роста
    };

    struct Edge {
        // стркутура ребра
        size_t from; // ID исходящей вершины
        size_t to; // ID входящей вершины
        double current; // ток ребра
        double sigma; // проводимость ребра
    };
    std::array<double, 3> start_r; // Координаты начала расчётной области
    std::array<double, 3> end_r; // Координаты конца расчётной области
    double h; // шаг сетки
    double delta_t; // шаг времени
    double r; // радиус канала
    double R; // радиу чехла
    double q_plus_max; // максимальный положительный заряд в точке
    double q_minus_max; // максимальный отрицательный заряд в точке
    double Q_plus_s; // максимальный положительный заряд в чехле
    double Q_minus_s; // максимальный отрицательный заряд в чехле
    double resistance; // Сопротивение
    double E_plus; // Положительная напряженность
    double E_minus; // ОТрицательная напряженность
    double alpha; // Параметр нагрезва звеньев
    double beta; // Параметр охлаждения звеньев
    double sigma; // Проводимость
    double degree_probability_growth; // Степень вероятности роста
    size_t periphery_size; // Максимальное количество переферийных зарядов
    std::function<double(const std::array<double, 3>&)> external_field_potential; // Функция внешнего поля напряжённости

    std::vector<Vertex> vertices; // вершины
    std::vector<bool> vertices_peripherality; // Переферийсность вершин
    std::vector<bool> vertices_activity; // Активность вершин
    std::vector<Edge> edges; // рёбра
    std::vector<bool> edges_activity; // активность рёбер

    DynamicOctree dynoctree; // структура хранения зарядов в виде восмиричного дерева

    // здесь индексы в массиве ребер

    std::vector<cubic_grid> graph; // граф дерева

    // struct Less {
    //     bool operator()(const std::array<int, 3>& lhs, const std::array<int, 3>& rhs) const {
    //         if (lhs[0] == rhs[0] && lhs[1] == rhs[1]) {
    //             return lhs[2] < rhs[2];
    //         }
    //         if (lhs[0] == rhs[0]) {
    //             return lhs[1] < rhs[1];
    //         }
    //         return lhs[0] < rhs[0];
    //     }
    // };

    // std::map<std::array<int, 3>, int, Less> internal_coords_to_id; // массив координат
    std::map<size_t, size_t> hash_to_id;
    size_t iter_number; // количество итераций

    int seed;
    int max_number_edges;
    mutable std::mt19937 gen;
    mutable std::uniform_real_distribution<> dis;
};
