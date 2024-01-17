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

        bool node_empty;
        double sum_q;
        std::array<double, 3> center_mass;

        std::array<double, 3> center;
        double size;
        
        std::vector<Node*> children;

        void add_charge(double charge, const std::array<double,3>& coords){
            //std::cout<<"Добавление заряда в октодерево\n";
            if(node_empty){
                sum_q = charge;
                center_mass = {coords[0], coords[1], coords[2]};
                node_empty = false;
            } else {
                if(children.empty()){
                    create_children();
                    size_t id_child = find_id_child(coords);
                    (*children[id_child]).add_charge(charge, coords);
                    
                    id_child = find_id_child(center_mass);
                    (*children[id_child]).add_charge(sum_q, center_mass);
                    
                    sum_q += charge;
                    // calc_sumCharge();
                    calc_center_mass();   
                } else {
                    size_t id_child = find_id_child(coords);
                    (*children[id_child]).add_charge(charge, coords);
                    sum_q += charge;
                    // calc_sumCharge();
                    calc_center_mass();                    
                }
            }

            // -- old --
            // if(children.empty()){
            //     if(node_empty){
            //         sum_q = charge;
            //         center_mass = {coords[0], coords[1], coords[2]};
            //         node_empty = false;              
            //     } else {
            //         create_children();
            //         Node child = find_id_child(coords);
            //         std::cout<<"Заряд до: "<<child.sum_q<<std::endl;
            //         child.add_charge(charge, coords);
            //         std::cout<<"Заряд после: "<<child.sum_q<<std::endl;
            //         child = find_id_child(center_mass);
            //         child.add_charge(sum_q, center_mass);
            //         sum_q += charge;
            //         // calc_sumCharge();
            //         calc_center_mass();                    
            //     }
            // } else {
            //     Node child = find_id_child(coords);
            //     child.add_charge(charge, coords);
            //     sum_q += charge;
            //     // calc_sumCharge();
            //     calc_center_mass();
            // }

        }

        double potencial_in_point(const std::array<double, 3>& point, double h){
            //std::cout<<"Расчёт потенциала в точке\n";
            if (node_empty) return 0;

            double x = center_mass[0];
            double y = center_mass[1];
            double z = center_mass[2];
            double l = std::sqrt(std::pow(point[0] - x, 2) + std::pow(point[1] - y, 2) + std::pow(point[2] - z, 2));
            double mirror_l = std::sqrt(std::pow(point[0] - x, 2) + std::pow(point[1] - y, 2) + std::pow(point[2] + z, 2));
            double k = 1 / (4 * std::numbers::pi * epsilon_0);
            double phi = 0;

            if(children.empty()){
                if (l < kEps) {
                    phi = sum_q * k * (1 / h - 1 / mirror_l);
                } else {
                    phi = sum_q * k * (1 / l - 1 / mirror_l);
                }
            } else {
                if(size < l * tetta){
                    phi = sum_q * k * (1 / l - 1 / mirror_l);
                } else {
                    for(auto& child : children){
                        phi += (*child).potencial_in_point(point, h);
                    }
                }
            }
            return phi;
        }

        void create_children(){
            //std::cout<<"Создание потомков\n";
            if(children.empty() && size > 100){
                for(int z = 0; z < 2; z++){
                    for(int y = 0; y < 2; y++){
                        for(int x = 0; x < 2; x++){
                            double cx = center[0] + (2*x - 1) * size / 4;
                            double cy = center[1] + (2*y - 1) * size / 4;
                            double cz = center[2] + (2*z - 1) * size / 4;
                            children.push_back(new Node{.node_empty = true,
                                                        .sum_q = 0,
                                                        .center_mass = {cx, cy, cz},
                                                        .center = {cx, cy, cz},
                                                        .size = size / 2,
                                                        .children = {}
                                                        });
                        }
                    }
                }
            }
        }

        size_t find_id_child(const std::array<double, 3>& coords){
            std::vector<size_t> temp = {0, 1, 2, 3, 4, 5, 6, 7};
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
                return temp[1];
            } else {
                return temp[0];
            }
        }
        
        void calc_center_mass(){
            std::array<double, 3> moment = {0, 0, 0};
            for(auto& child : children){
                if(!(*child).node_empty){
                    moment[0] += (*child).center_mass[0] * (*child).sum_q;
                    moment[1] += (*child).center_mass[1] * (*child).sum_q;
                    moment[2] += (*child).center_mass[2] * (*child).sum_q;
                }
            }
            if(sum_q != 0){
                center_mass = {moment[0] / (sum_q),
                                moment[1] / (sum_q),
                                moment[2] / (sum_q)};
            } else {
                center_mass = {center[0], center[1], center[2]};
            }
        }

        double get_charge(const std::array<double, 3>& point){
            double charge;
            if (children.empty()){
                if(!node_empty) charge = sum_q;
                else std::cout<<"Ошибка поиска заряда\n";
            } else {
                size_t id = find_id_child(point);
                charge = (*children[id]).get_charge(point);
            }
            return charge;
        }

        friend std::ostream& operator<<(std::ostream& stream, const Node& node)
        {
            stream << "Node: "<< &node;
            stream << "\n\tnode_empty: "<<node.node_empty;
            stream << "\n\tsum_q: "<<node.sum_q;
            stream << "\n\tcenter_mass: ("<<node.center_mass[0]<<", "<<node.center_mass[1]<<", "<<node.center_mass[2]<<")";
            stream << "\n\tcenter: ("<<node.center[0]<<", "<<node.center[1]<<", "<<node.center[2]<<")";
            stream << "\n\tsize: "<<node.size;
            stream << "\n\tchildren: ";
            if (node.children.empty()){
                stream << "Нет потомков"<<"\n";
            } else {
                stream << "[";
                for(auto& child : node.children){
                    stream << child <<", ";
                }
                stream << "]\n";
            }
            return stream;
        }

        ~Node(){
            if(!children.empty()){
                for(auto& child : children){
                    delete child;
                }
            }
        }
        
    } root;

    public:
    Octree(){}

    Octree(const std::array<double, 3>& start, const std::array<double, 3>& end) {
        root = Node{.node_empty = true,
                    .sum_q = 0,
                    .center_mass = {(end[0] + start[0]) / 2, (end[1] + start[1]) / 2, (end[2] + start[2]) / 2},
                    .center = {(end[0] + start[0]) / 2, (end[1] + start[1]) / 2, (end[2] + start[2]) / 2},
                    .size = std::max(std::max(end[0] - start[0], end[1] - start[1]), std::max(end[1] - start[1], end[2] - start[2])),
                    .children = {}};
    }

    void add_charge(double charge, const std::array<double, 3>& coords) {
        // std::cout<<"Добавление заряда в октодерево\n";
        root.add_charge(charge, coords);
    }

    double potencial_in_point(const std::array<double, 3>& point, double h) {
        //std::cout<<"Расчёт потенциала в точке\n";
        return root.potencial_in_point(point, h);
    }

    double get_charge(const std::array<double, 3>& point){
        return root.get_charge(point);
    }

    friend std::ostream& operator<<(std::ostream& stream, const Octree& octree)
    {
        stream << octree.root <<"\n";
        return stream;
    }

    ~Octree(){}
};