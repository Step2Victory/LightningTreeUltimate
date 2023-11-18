#pragma once
#include <iostream>
#include <vector>
#include <array>
#include <optional>
#include <memory>
#include <cmath>
#include <numbers>
#include "Constants.h"

constexpr double tetta = 0.1;

class Octree{

    private:

    struct Node {

        std::optional<size_t> id_charge;
        double sum_q;
        double sum_Q;
        std::optional<std::array<double, 3>> center_mass;

        std::array<double, 3> center;
        double size;
        
        std::vector<std::unique_ptr<Node>> children;

        void add_charge(size_t id, const auto& charges){
            //std::cout<<"Добавление заряда в октодерево\n";
            if(children.empty()){
                if(id_charge.has_value()){
                    create_children();
                    (find_child(id, charges)).add_charge(id, charges);
                    (find_child(*id_charge, charges)).add_charge(*id_charge, charges);
                    id_charge = NULL;
                    sum_q += charges[id].q;
                    sum_Q += charges[id].Q;
                    // calc_sumCharge();
                    calc_center_mass();
                } else {
                    sum_q = charges[id].q;
                    sum_Q = charges[id].Q;
                    center_mass = {charges[id].coords[0], charges[id].coords[1], charges[id].coords[2]};
                }
            } else {
                (find_child(id, charges)).add_charge(id, charges);
                sum_q += charges[id].q;
                sum_Q += charges[id].Q;
                // calc_sumCharge();
                calc_center_mass();
            }
        }

        void add_charges(auto& charges){
            for(size_t i = 0; i < charges.size(); i++){
                add_charge(i, charges);
            }
        }

        void delete_charge(size_t id, const auto& charges){
            //std::cout<<"Удаление заряда из октодерева\n";
            if(id_charge){
                if(id_charge == id) this->clear();
                else std::cout<<"Ошибка поиска!\n";
            } else {
                (find_child(id, charges)).delete_charge(id, charges);
                if(!delete_children()){
                    this->sum_q -= charges[id].q;
                    this->sum_Q -= charges[id].Q;
                }
                this->calc_center_mass();
            }
            
        }

        void recalc_sumCharge(const auto& charges){
            //std::cout<<"Перерасчет суммы зарядов\n";
            if(children.empty()){
                if(id_charge.has_value()){
                    this->sum_q = charges[*id_charge].q;
                    this->sum_Q = charges[*id_charge].Q;
                }
            } else {
                double q = 0, Q = 0;
                std::array<double, 3> moment = {0, 0, 0};
                for(auto& child : children){
                    (*child).recalc_sumCharge(charges);
                    q += (*child).sum_q;
                    Q += (*child).sum_Q;
                    moment[0] += ((*child).sum_q + (*child).sum_Q) * (*child).center_mass.value()[0];
                    moment[1] += ((*child).sum_q + (*child).sum_Q) * (*child).center_mass.value()[1];
                    moment[2] += ((*child).sum_q + (*child).sum_Q) * (*child).center_mass.value()[2];
                }
                sum_q = q;
                sum_Q = Q;
                center_mass = {moment[0] / (q + Q), moment[1] / (q + Q), moment[2] / (q + Q)};
            }
        }

        double potencial_in_point(const std::array<double, 3>& point, double r, double R, double h){
            //std::cout<<"Расчёт потенциала в точке\n";
            if (!center_mass.has_value()) return 0;

            double x = (*center_mass)[0];
            double y = (*center_mass)[1];
            double z = (*center_mass)[2];
            double l = std::sqrt(std::pow(point[0] - x, 2) + std::pow(point[1] - y, 2) + std::pow(point[2] - z, 2));
            double mirror_l = std::sqrt(std::pow(point[0] - x, 2) + std::pow(point[1] - y, 2) + std::pow(point[2] + z, 2));
            double k = 1 / (4 * std::numbers::pi * epsilon_0);

            if(id_charge.has_value()){
                if (l < kEps) {
                    return sum_q * k * (1 / (h / 2 + r) - 1 / (mirror_l + r)) +
                        sum_Q * k * (1 / (h / 2 + R) - 1 / (mirror_l + R));
                } else {
                    return sum_q * k * (1 / (l + r) - 1 / (mirror_l + r)) +
                        sum_Q * k * (1 / (l + R) - 1 / (mirror_l + R));
                }
            }

            if(size / l < tetta){
                if (l < kEps) {
                    return sum_q * k * (1 / (h / 2 + r) - 1 / (mirror_l + r)) +
                        sum_Q * k * (1 / (h / 2 + R) - 1 / (mirror_l + R));
                } else {
                    return sum_q * k * (1 / (l + r) - 1 / (mirror_l + r)) +
                        sum_Q * k * (1 / (l + R) - 1 / (mirror_l + R));
                }
            } else {
                double phi = 0;
                for(auto& child : children){
                    phi += (*child).potencial_in_point(point, r, R, h);
                }
                return phi;
            }
        }

        void create_children(){
            //std::cout<<"Создание потомков\n";
            if(children.empty()){
                for(int i = 0; i < 2; i++){
                    for(int j = 0; j < 2; j++){
                        for(int k = 0; k < 2; k++){
                            children.push_back(std::make_unique<Node>(
                                        Node{.id_charge = NULL,
                                             .sum_q = 0,
                                             .sum_Q = 0,
                                             .center = {center[0] + (2*i - 1) * size / 4,
                                                        center[1] + (2*j - 1) * size / 4,
                                                        center[2] + (2*k - 1) * size / 4},
                                             .size = size / 2}));
                        }
                    }
                }
            }
        }

        bool delete_children(){
            //std::cout<<"Удаление потомков\n";
            int count = 0;
            size_t id;
            for(size_t i = 0; i < children.size(); i++){
                if((*children[i]).center_mass.has_value()) {
                    count++;
                    id = i;
                }
            }
            if(count == 1){
                this->id_charge = (*children[id]).id_charge;
                this->sum_q = (*children[id]).sum_q;
                this->sum_Q = (*children[id]).sum_Q;
                this->center_mass = (*children[id]).center_mass;
                children.clear();
                return true;
            }
            if(count == 0) {
                children.clear();
                return true;
            }
            return false;
        }

        void clear(){
            id_charge = NULL;
            sum_q = 0;
            sum_Q = 0;
            center_mass = {NULL, NULL, NULL};
        }

        Node find_child(size_t id, const auto& charges){
            std::vector<std::unique_ptr<Node>> temp(children);
            if(charges[id].coords[2] > center[2]){
                temp.erase(temp.begin(), temp.begin() + 4);
            } else {
                temp.erase(temp.begin() + 4, temp.end());
            }

            if(charges[id].coords[1] > center[1]){
                temp.erase(temp.begin(), temp.begin() + 2);
            } else {
                temp.erase(temp.begin() + 2, temp.end());
            }

            if(charges[id].coords[0] > center[0]){
                return *temp[1];
            } else {
                return *temp[0];
            }
        }
        
        // void calc_sumCharge(){
        //     double q = 0;
        //     double Q = 0;
        //     for(auto& child : children){
        //         q += child.sum_q;
        //         Q += child.sum_Q;
        //     } 
        //     this->sum_q = q;
        //     this->sum_Q = Q;
        // }

        void calc_center_mass(){
            std::array<double, 3> moment = {0, 0, 0};
            for(auto& child : children){
                if(double sum_charge = ((*child).sum_Q + (*child).sum_q) != 0){
                    moment[0] += (*child).center_mass.value()[0] * (sum_charge);
                    moment[1] += (*child).center_mass.value()[1] * (sum_charge);
                    moment[2] += (*child).center_mass.value()[2] * (sum_charge);
                }
            }
            if(double sum_charge = (sum_Q + sum_q) != 0){
                center_mass = {moment[0] / (sum_charge),
                                moment[1] / (sum_charge),
                                moment[2] / (sum_charge)};
            }
        }
    } root;

    public:
    Octree(){}

    Octree(const std::array<double, 3>& start, const std::array<double, 3>& end) {
        root = Node{.id_charge = NULL,
                    .sum_q = 0,
                    .sum_Q = 0,
                    .center = {(end[0] + start[0]) / 2, (end[1] + start[1]) / 2, (end[2] + start[2]) / 2},
                    .size = std::max(std::max(end[0] - start[0], end[1] - start[1]), std::max(end[1] - start[1], end[2] - start[2]))};
    }

    Octree(const std::array<double, 3>& start, const std::array<double, 3>& end, const auto& charges) {
        root = Node{.id_charge = NULL,
                    .sum_q = 0,
                    .sum_Q = 0,
                    .center = {(end[0] + start[0]) / 2, (end[1] + start[1]) / 2, (end[2] + start[2]) / 2},
                    .size = std::max(std::max(end[0] - start[0], end[1] - start[1]), std::max(end[1] - start[1], end[2] - start[2]))};
        root.add_charges(charges);
    }

    void add_charge(size_t id, const auto& charges) {
        // std::cout<<"Добавление заряда в октодерево\n";
        root.add_charge(id, charges);
    }

    void delete_charge(size_t id, const auto& charges){
        // std::cout<<"Удаление заряда из октодерева\n";
        root.delete_charge(id, charges);
    }

    void recalc_sumCharge(const auto& charges) {
        // std::cout<<"Перерасчет суммы зарядов\n";
        root.recalc_sumCharge(charges);
    }

    double potencial_in_point(const std::array<double, 3>& point, double r, double R, double h) {
        //std::cout<<"Расчёт потенциала в точке\n";
        return root.potencial_in_point(point, r, R, h);
    }

    ~Octree(){}
};