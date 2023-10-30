#pragma once
#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <numbers>
#include "Constants.h"

constexpr double tetta = 0.1;

class Octree{

    private:

    struct Node {

        size_t id_charge;
        double sum_q;
        double sum_Q;
        std::array<double, 3> center_mass;

        std::array<double, 3> center;
        double size;
        
        std::vector<Node> children;

        void add_charge(size_t id, const auto& charges){
            if(children.empty()){
                if(id_charge){
                    create_children();
                    (find_child(id, charges)).add_charge(id, charges);
                    (find_child(id_charge, charges)).add_charge(id_charge, charges);
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
            if(children.empty()){
                if(id_charge){
                    this->sum_q = charges[id_charge].q;
                    this->sum_Q = charges[id_charge].Q;
                }
            } else {
                double q = 0, Q = 0;
                std::array<double, 3> moment = {0, 0, 0};
                for(auto& child : children){
                    child.recalc_sumCharge(charges);
                    q += child.sum_q;
                    Q += child.sum_Q;
                    moment[0] += (child.sum_q + child.sum_Q) * child.center_mass[0];
                    moment[1] += (child.sum_q + child.sum_Q) * child.center_mass[1];
                    moment[2] += (child.sum_q + child.sum_Q) * child.center_mass[2];
                }
                this->sum_q = q;
                this->sum_Q = Q;
                this->center_mass = {moment[0] / (q + Q), moment[1] / (q + Q), moment[2] / (q + Q)};
            }
        }

        double potencial_in_point(const std::array<double, 3>& point, double r, double R, double h){
            if (center_mass.empty()) return 0;

            double x = this->center_mass[0];
            double y = this->center_mass[1];
            double z = this->center_mass[2];
            double l = std::sqrt(std::pow(point[0] - x, 2) + std::pow(point[1] - y, 2) + std::pow(point[2] - z, 2));
            double mirror_l = std::sqrt(std::pow(point[0] - x, 2) + std::pow(point[1] - y, 2) + std::pow(point[2] + z, 2));
            double k = 1 / (4 * std::numbers::pi * epsilon_0);

            if(this->id_charge){
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
                    phi += child.potencial_in_point(point, r, R, h);
                }
                return phi;
            }
        }

        void create_children(){
            if(children.empty()){
                for(int i = -1; i < 2; i+=2){
                    for(int j = -1; j < 2; j+=2){
                        for(int k = -1; k < 2; k+=2){
                            children.push_back(Node{.id_charge = NULL,
                                                    .sum_q = 0,
                                                    .sum_Q = 0,
                                                    .center = {center[0] + i * size / 4,
                                                                center[1] + j * size / 4,
                                                                center[2] + k * size / 4},
                                                    .size = size / 2});
                        }
                    }
                }
            }
        }

        bool delete_children(){
            int count = 0;
            size_t id;
            for(size_t i = 0; i < children.size(); i++){
                if(!children[i].center_mass.empty()) {
                    count++;
                    id = i;
                }
            }
            if(count == 1){
                this->id_charge = children[id].id_charge;
                this->sum_q = children[id].sum_q;
                this->sum_Q = children[id].sum_Q;
                this->center_mass = children[id].center_mass;
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
            this->center_mass = std::array<double, 3>{};
        }

        Node find_child(size_t id_charge, const auto& charges){
            std::vector<Node> temp(children);
            if(charges[id_charge].coords[2] > center[2]){
                temp.erase(temp.begin(), temp.begin() + 4);
            } else {
                temp.erase(temp.begin() + 4, temp.end());
            }

            if(charges[id_charge].coords[1] > center[1]){
                temp.erase(temp.begin(), temp.begin() + 2);
            } else {
                temp.erase(temp.begin() + 2, temp.end());
            }

            if(charges[id_charge].coords[0] > center[0]){
                return temp[1];
            } else {
                return temp[0];
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
                if(double sum_charge = (child.sum_Q + child.sum_q) != 0){
                    moment[0] += child.center_mass[0] * (sum_charge);
                    moment[1] += child.center_mass[1] * (sum_charge);
                    moment[2] += child.center_mass[2] * (sum_charge);
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
        root.add_charge(id, charges);
    }

    void delete_charge(size_t id, const auto& charges){
        root.delete_charge(id, charges);
    }

    void recalc_sumCharge(const auto& charges) {
        root.recalc_sumCharge(charges);
    }

    double potencial_in_point(const std::array<double, 3>& point, double r, double R, double h) {
        return root.potencial_in_point(point, r, R, h);
    }

    ~Octree(){}
};