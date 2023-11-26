#pragma once
#include <iostream>
#include <vector>
#include <array>
#include <optional>
#include <memory>
#include <cmath>
#include <numbers>
#include "Constants.h"

// constexpr double tetta = 0.1;

class Octree{

    private:

    struct Node {

        double sum_q;
        std::optional<std::array<double, 3>> center_mass;

        std::array<double, 3> center;
        double size;
        
        std::vector<Node*> children;

        void add_charge(double charge, const std::array<double,3>& coords){
            //std::cout<<"Добавление заряда в октодерево\n";
            if(children.empty()){
                if(center_mass.has_value()){
                    create_children();
                    (find_child(coords)).add_charge(charge, coords);
                    (find_child(*center_mass)).add_charge(sum_q, *center_mass);
                    sum_q += charge;
                    // calc_sumCharge();
                    calc_center_mass();
                } else {
                    sum_q = charge;
                    center_mass = {coords[0], coords[1], coords[2]};
                }
            } else {
                (find_child(coords)).add_charge(charge, coords);
                sum_q += charge;
                // calc_sumCharge();
                calc_center_mass();
            }
        }

        double potencial_in_point(const std::array<double, 3>& point, double h){
            //std::cout<<"Расчёт потенциала в точке\n";
            if (!center_mass.has_value()) return 0;

            double x = (*center_mass)[0];
            double y = (*center_mass)[1];
            double z = (*center_mass)[2];
            double l = std::sqrt(std::pow(point[0] - x, 2) + std::pow(point[1] - y, 2) + std::pow(point[2] - z, 2));
            double mirror_l = std::sqrt(std::pow(point[0] - x, 2) + std::pow(point[1] - y, 2) + std::pow(point[2] + z, 2));
            double k = 1 / (4 * std::numbers::pi * epsilon_0);

            if(children.empty()){
                if (l < kEps) {
                    return sum_q * k * (1 / h - 1 / mirror_l);
                } else {
                    return sum_q * k * (1 / l - 1 / mirror_l);
                }
            }
            double distance = std::sqrt(std::pow(point[0] - center[0], 2) + std::pow(point[1] - center[1], 2) + std::pow(point[2] - center[2], 2));
            if(size / distance < tetta){
                if (l < kEps) {
                    return sum_q * k * (1 / h - 1 / mirror_l);
                } else {
                    return sum_q * k * (1 / l - 1 / mirror_l);
                }
            } else {
                double phi = 0;
                for(auto& child : children){
                    phi += (*child).potencial_in_point(point, h);
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
                            children.push_back(new Node{.sum_q = 0,
                                             .center = {center[0] + (2*i - 1) * size / 4,
                                                        center[1] + (2*j - 1) * size / 4,
                                                        center[2] + (2*k - 1) * size / 4},
                                             .size = size / 2});
                        }
                    }
                }
            }
        }

        Node find_child(const std::array<double, 3>& coords){
            std::vector<Node*> temp(children);
            if(coords[2] > center[2]){
                temp.erase(temp.begin(), temp.begin() + 4);
            } else {
                temp.erase(temp.begin() + 4, temp.end());
            }

            if(coords[1] > center[1]){
                temp.erase(temp.begin(), temp.begin() + 2);
            } else {
                temp.erase(temp.begin() + 2, temp.end());
            }

            if(coords[0] > center[0]){
                return *temp[1];
            } else {
                return *temp[0];
            }
        }
        
        void calc_center_mass(){
            std::array<double, 3> moment = {0, 0, 0};
            for(auto& child : children){
                if((*child).sum_q != 0){
                    moment[0] += (*child).center_mass.value()[0] * (sum_q);
                    moment[1] += (*child).center_mass.value()[1] * (sum_q);
                    moment[2] += (*child).center_mass.value()[2] * (sum_q);
                }
            }
            if(sum_q != 0){
                center_mass = {moment[0] / (sum_q),
                                moment[1] / (sum_q),
                                moment[2] / (sum_q)};
            }
        }

        ~Node(){
            for(auto& child : children){
                delete child;
            }
        }
        
    } root;

    public:
    Octree(){}

    Octree(const std::array<double, 3>& start, const std::array<double, 3>& end) {
        root = Node{.sum_q = 0,
                    .center = {(end[0] + start[0]) / 2, (end[1] + start[1]) / 2, (end[2] + start[2]) / 2},
                    .size = std::max(std::max(end[0] - start[0], end[1] - start[1]), std::max(end[1] - start[1], end[2] - start[2]))};
    }

    void add_charge(double charge, const std::array<double, 3>& coords) {
        // std::cout<<"Добавление заряда в октодерево\n";
        root.add_charge(charge, coords);
    }

    double potencial_in_point(const std::array<double, 3>& point, double h) {
        //std::cout<<"Расчёт потенциала в точке\n";
        return root.potencial_in_point(point, h);
    }

    ~Octree(){
        root.~Node();
    }
};