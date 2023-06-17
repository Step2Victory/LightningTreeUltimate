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
    // TO DO
}

void LightningTree::CountElectricity() {
    // TO DO
}

void LightningTree::CountCurrent() {
    // TO DO
}

void LightningTree::Transport() {
    // TO DO
}

void LightningTree::Grow() {
    // TO DO
}

void LightningTree::Delete() {
    // TO DO
}

bool LightningTree::GrowthCriterion(
    size_t vertex_id, const std::array<double, 3>& coords) const {
    // TO DO
}

bool LightningTree::DeletionCriterion(size_t vertex_id) const {
    // TO DO
}