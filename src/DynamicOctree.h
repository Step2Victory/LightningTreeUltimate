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

class DynamicOctree{

    private:

    struct DynamicNode {

        bool node_empty;
        size_t id_charge;
        double sum_q;
        double sum_Q;
        std::array<double, 3> center_mass;

        std::array<double, 3> center;
        double size;
        int lvl;
        
        std::vector<std::unique_ptr<DynamicNode>> children;

        void add_charge(size_t id, const auto& charges){
            // std::cout<<"Добавление заряда id: "<<id<<", координаты {"<<charges[id].coords[0]<<", "<<charges[id].coords[1]<<", "<<charges[id].coords[2]<<"}"<<" в октодерево, уровень: "<<lvl<<'\n';
            if(node_empty){
                // std::cout<<"Добавление в пустую ноду"<<'\n';
                node_empty = false;
                id_charge = id;
                sum_q = charges[id].q;
                sum_Q = charges[id].Q;
                center_mass = {charges[id].coords[0], charges[id].coords[1], charges[id].coords[2]};
            } else {
                if(children.empty() && create_children()){
                    // std::cout<<"Добавление в новую ноду\n";
                    size_t id_child = find_id_child(charges[id].coords);
                    // std::cout<<"id потомка: "<<id_child<<" с координатами центра {"<<(*children[id_child]).center[0]<<", "<<(*children[id_child]).center[1]<<", "<<(*children[id_child]).center[2]<<"} и размером ноды: "<<(*children[id_child]).size<<'\n';
                    (*children[id_child]).add_charge(id, charges);
                    // std::cout<<"добавление заряда в ноду: "<<children[id_child]<<" с id заряда: "<<children[id_child]->id_charge<<", переданный id "<<id<<'\n';
                    id_child = find_id_child(charges[id_charge].coords);
                    // std::cout<<"id потомка: "<<id_child<<" с координатами центра {"<<(*children[id_child]).center[0]<<", "<<(*children[id_child]).center[1]<<", "<<(*children[id_child]).center[2]<<"} и размером ноды: "<<(*children[id_child]).size<<'\n';
                    (*children[id_child]).add_charge(id_charge, charges);
                    // std::cout<<"добавление заряда в ноду: "<<children[id_child]<<" с id заряда: "<<children[id_child]->id_charge<<", переданный id "<<id_charge<<'\n';
                    id_charge = 0;
                    sum_q += charges[id].q;
                    sum_Q += charges[id].Q;
                    // std::cout<<"расчёт центра масс после создания потомков\n";
                    calc_center_mass();
                    // std::cout<<"расчёт центра масс завершился\n";
                } else {
                    // std::cout<<"Поиск свободной ноды и добавление заряда\n";
                    size_t id_child = find_id_child(charges[id].coords);
                    // std::cout<<"id потомка: "<<id_child<<" с координатами центра {"<<(*children[id_child]).center[0]<<", "<<(*children[id_child]).center[1]<<", "<<(*children[id_child]).center[2]<<"} и размером ноды: "<<(*children[id_child]).size<<'\n';
                    (*children[id_child]).add_charge(id, charges);
                    sum_q += charges[id].q;
                    sum_Q += charges[id].Q;
                    // std::cout<<"расчёт центра масс после добавления заряда\n";
                    calc_center_mass();
                    // std::cout<<"расчёт центра масс завершился\n";            
                }
                // print('d');
            }

            // -- old --
            // if(children.empty()){
            //     if(id_charge.has_value()){
            //         create_children();
            //         (find_id_child(id, charges)).add_charge(id, charges);
            //         (find_id_child(*id_charge, charges)).add_charge(*id_charge, charges);
            //         id_charge = NULL;
            //         sum_q += charges[id].q;
            //         sum_Q += charges[id].Q;
            //         // calc_sumCharge();
            //         calc_center_mass();
            //     } else {
            //         sum_q = charges[id].q;
            //         sum_Q = charges[id].Q;
            //         center_mass = {charges[id].coords[0], charges[id].coords[1], charges[id].coords[2]};
            //     }
            // } else {
            //     (find_id_child(id, charges)).add_charge(id, charges);
            //     sum_q += charges[id].q;
            //     sum_Q += charges[id].Q;
            //     // calc_sumCharge();
            //     calc_center_mass();
            // }
        }

        // -- old --
        // void add_charges(auto& charges){
        //     for(size_t i = 0; i < charges.size(); i++){
        //         add_charge(i, charges);
        //     }
        // }

        void delete_charge(size_t id, const auto& charges){
            //std::cout<<"Удаление заряда из октодерева\n";
            if(children.empty()){
                if(node_empty){
                    std::cout<<"Ошибка поиска! Пустая нода\n";
                    // print('d');
                } else {
                    if(id_charge == id) {
                        // print('d');
                        clear();
                    }
                    else std::cout<<"Ошибка поиска! Не совпадение id вершин\n";
                }
            } else {
                size_t id_child = find_id_child(charges[id].coords);
                (*children[id_child]).delete_charge(id, charges);
                if(!delete_children()){
                    sum_q -= charges[id].q;
                    sum_Q -= charges[id].Q;
                    // std::cout<<"расчёт центра масс после удаления потомков\n";
                    calc_center_mass();
                    // std::cout<<"расчёт центра масс завершился\n";
                    // print('d');
                }
            }
            // -- old --
            // if(!node_empty){
            //     if(id_charge == id) clear();
            //     else std::cout<<"Ошибка поиска!\n";
            // } else {
            //     (find_id_child(id, charges)).delete_charge(id, charges);
            //     if(!delete_children()){
            //         sum_q -= charges[id].q;
            //         sum_Q -= charges[id].Q;
            //         calc_center_mass();
            //     }
            // }
            
        }

        void recalc_sumCharge(const auto& charges){
            //std::cout<<"Перерасчет суммы зарядов\n";
            if(children.empty()){
                if(!node_empty){
                    sum_q = charges[id_charge].q;
                    sum_Q = charges[id_charge].Q;
                }
            } else {
                double q = 0, Q = 0;
                std::array<double, 3> moment = {0, 0, 0};
                for(auto& child : children){
                    (*child).recalc_sumCharge(charges);
                    q += (*child).sum_q;
                    Q += (*child).sum_Q;
                    moment[0] += ((*child).sum_q + (*child).sum_Q) * (*child).center_mass[0];
                    moment[1] += ((*child).sum_q + (*child).sum_Q) * (*child).center_mass[1];
                    moment[2] += ((*child).sum_q + (*child).sum_Q) * (*child).center_mass[2];
                }
                sum_q = q;
                sum_Q = Q;
                double sum_charge = (sum_Q + sum_q);
                if (sum_charge > kEps || sum_charge < -kEps)
                    center_mass = {moment[0] / sum_charge, moment[1] / sum_charge, moment[2] / sum_charge};
                else 
                    center_mass = center;
            }
            // std::cout<<"Результат перерасчёта\n";
            // print('p');
        }

        double potencial_in_point(const std::array<double, 3>& point, double r, double R, double h){
            //std::cout<<"Расчёт потенциала в точке\n";
            if (node_empty) {
                // std::cout<<"Пустая нода\n";
                return 0;
            }

            double x = center_mass[0];
            double y = center_mass[1];
            double z = center_mass[2];
            double l = std::sqrt(std::pow(point[0] - x, 2) + std::pow(point[1] - y, 2) + std::pow(point[2] - z, 2));
            double mirror_l = std::sqrt(std::pow(point[0] - x, 2) + std::pow(point[1] - y, 2) + std::pow(point[2] + z, 2));
            double k = 1 / (4 * std::numbers::pi * epsilon_0);
            double phi = 0;

            if(children.empty()){
                if (l < kEps) {
                    phi = sum_q * k * (1 / (h / 2 + r) - 1 / (mirror_l + r)) + 
                            sum_Q * k * (1 / (h / 2 + R) - 1 / (mirror_l + R));
                } else {
                    phi = sum_q * k * (1 / (l + r) - 1 / (mirror_l + r)) +
                            sum_Q * k * (1 / (l + R) - 1 / (mirror_l + R));
                }
            } else {
                if(size < l * tetta){
                    phi = sum_q * k * (1 / (l + r) - 1 / (mirror_l + r)) +
                            sum_Q * k * (1 / (l + R) - 1 / (mirror_l + R));
                } else {
                    for(auto& child : children){
                        phi += (*child).potencial_in_point(point, r, R, h);
                    }
                }
            }
            return phi;

            // -- old --
            // if(id_charge.has_value()){
            //     if (l < kEps) {
            //         return sum_q * k * (1 / (h / 2 + r) - 1 / (mirror_l + r)) +
            //             sum_Q * k * (1 / (h / 2 + R) - 1 / (mirror_l + R));
            //     } else {
            //         return sum_q * k * (1 / (l + r) - 1 / (mirror_l + r)) +
            //             sum_Q * k * (1 / (l + R) - 1 / (mirror_l + R));
            //     }
            // }
            // if(size / l < tetta){
            //     if (l < kEps) {
            //         return sum_q * k * (1 / (h / 2 + r) - 1 / (mirror_l + r)) +
            //             sum_Q * k * (1 / (h / 2 + R) - 1 / (mirror_l + R));
            //     } else {
            //         return sum_q * k * (1 / (l + r) - 1 / (mirror_l + r)) +
            //             sum_Q * k * (1 / (l + R) - 1 / (mirror_l + R));
            //     }
            // } else {
            //     double phi = 0;
            //     for(auto& child : children){
            //         phi += (*child).potencial_in_point(point, r, R, h);
            //     }
            //     return phi;
            // }
        }

        bool create_children(){
            // std::cout<<"Создание потомков\n";
            if(children.empty() && size > 100){
                for(int z = 0; z < 2; z++){
                    for(int y = 0; y < 2; y++){
                        for(int x = 0; x < 2; x++){
                            double cx = center[0] + (2*x - 1) * size / 4;
                            double cy = center[1] + (2*y - 1) * size / 4;
                            double cz = center[2] + (2*z - 1) * size / 4;
                            children.push_back(std::make_unique<DynamicNode>(
                                                    DynamicNode{.node_empty = true,
                                                                .id_charge = 0,
                                                                .sum_q = 0,
                                                                .sum_Q = 0,
                                                                .center_mass = {cx, cy, cz},
                                                                .center = {cx, cy, cz},
                                                                .size = size / 2,
                                                                .lvl = lvl + 1,
                                                                .children = {}
                                                                }));
                        }
                    }
                }
            return true;
            }
            return false;
        }

        bool delete_children(){
            //std::cout<<"Удаление потомков\n";
            int count = 0;
            size_t id;
            for(size_t i = 0; i < children.size(); i++){
                if(!node_empty) {
                    count++;
                    id = i;
                }
            }
            if(count == 1){
                id_charge = (*children[id]).id_charge;
                sum_q = (*children[id]).sum_q;
                sum_Q = (*children[id]).sum_Q;
                center_mass = (*children[id]).center_mass;
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
            node_empty = true;
            id_charge = 0;
            sum_q = 0;
            sum_Q = 0;
            center_mass = center;
            children.clear();
        }

        size_t find_id_child(const auto& coords){
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
        
        // void calc_sumCharge(){
        //     double q = 0;
        //     double Q = 0;
        //     for(auto& child : children){
        //         q += child.sum_q;
        //         Q += child.sum_Q;
        //     } 
        //     sum_q = q;
        //     sum_Q = Q;
        // }

        void calc_center_mass(){
            // std::cout<<"Запуск расчёта центра масс\n";
            std::array<double, 3> moment = {0, 0, 0};
            for(auto& child : children){
                if(!(*child).node_empty){
                    moment[0] += (*child).center_mass[0] * ((*child).sum_Q + (*child).sum_q);
                    moment[1] += (*child).center_mass[1] * ((*child).sum_Q + (*child).sum_q);
                    moment[2] += (*child).center_mass[2] * ((*child).sum_Q + (*child).sum_q);
                }
            }
            double sum_charge = (sum_Q + sum_q);
            // std::cout<<"сумма зарядов = "<<sum_charge<<'\n';
            if(sum_charge > kEps || sum_charge < -kEps){
                center_mass = {moment[0] / (sum_charge),
                                moment[1] / (sum_charge),
                                moment[2] / (sum_charge)};
            } else {
                center_mass = center;
            }
        }

        void print(char config)
        {
            std::cout << "DynamicNode: "<< this;
            std::cout << "\n\tnode_empty: "<<node_empty;
            std::cout << "\n\tid_charge: "<<id_charge;
            std::cout << "\n\tsum_q: "<<sum_q;
            std::cout << "\n\tsum_Q: "<<sum_Q;
            std::cout << "\n\tcenter_mass: ("<<center_mass[0]<<", "<<center_mass[1]<<", "<<center_mass[2]<<")";
            std::cout << "\n\tcenter: ("<<center[0]<<", "<<center[1]<<", "<<center[2]<<")";
            std::cout << "\n\tsize: "<<size;
            std::cout << "\n\tchildren: ";
            if (children.empty()){
                std::cout << "Нет потомков"<<"\n";
            } else {
                std::cout << "[";
                if(config == 'd') {
                    for(auto& child : children){
                        (*child).print(config);
                        std::cout << ", ";
                    }
                } else {
                    for(auto& child : children){
                        std::cout << &child << ", ";
                    }
                }
                std::cout << "]\n";
            }
        }

    } root;

    public:
    DynamicOctree(){}

    DynamicOctree(const std::array<double, 3>& start, const std::array<double, 3>& end) {
        root = DynamicNode{.node_empty = true,
                            .id_charge = 0,
                            .sum_q = 0,
                            .sum_Q = 0,
                            .center_mass = {(end[0] + start[0]) / 2, (end[1] + start[1]) / 2, (end[2] + start[2]) / 2},
                            .center = {(end[0] + start[0]) / 2, (end[1] + start[1]) / 2, (end[2] + start[2]) / 2},
                            .size = std::max(std::max(end[0] - start[0], end[1] - start[1]), std::max(end[1] - start[1], end[2] - start[2])),
                            .lvl = 0,
                            .children = {}
                            };
    }

    // -- old --
    // DynamicOctree(const std::array<double, 3>& start, const std::array<double, 3>& end, const auto& charges) {
    //     root = DynamicNode{.node_empty = true,
    //                         .id_charge = 0,
    //                         .sum_q = 0,
    //                         .sum_Q = 0,
    //                         .center_mass = {(end[0] + start[0]) / 2, (end[1] + start[1]) / 2, (end[2] + start[2]) / 2},
    //                         .center = {(end[0] + start[0]) / 2, (end[1] + start[1]) / 2, (end[2] + start[2]) / 2},
    //                         .size = std::max(std::max(end[0] - start[0], end[1] - start[1]), std::max(end[1] - start[1], end[2] - start[2])),
    //                         .children = {}
    //                         };
    //     root.add_charges(charges);
    // }

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

    void print()
    {
        root.print('p');
    }

    void print(char config)
    {
        root.print(config);
    }
};