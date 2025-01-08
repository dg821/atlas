#pragma once

#include <Eigen/Dense>
#include <vector>
#include <functional>
#include <string>
#include <iostream>
#include "NumericalIntegrator.h"  // Assuming this base class exists

class Dopri87 : public NumericalIntegrator {
public:
    using Vector6d = Eigen::Matrix<double, 6, 1>;
    using StateVector = Vector6d;
    using RightHandSide = std::function<StateVector(double, const StateVector&)>;

    Dopri87();
    ~Dopri87() override = default;

    StateVector step(const RightHandSide& f, double t, const StateVector& y, double h) const override;

    std::vector<std::pair<double, StateVector>> integrate(
        const RightHandSide& f, double t_start, double t_end, const StateVector& y_start,
        double h, double h_min, double h_max) const override;

    int getOrder() const override;
    std::string getName() const override;

private:
    static constexpr int stages = 13;
    static constexpr int order = 8;

    Eigen::Matrix<double, stages, stages> a;
    Eigen::Vector<double, stages> c;
    Eigen::Vector<double, stages> b7, b8;

    void initializeCoefficients();

    std::pair<StateVector, StateVector> computeStage(
        const RightHandSide& f, double t, const StateVector& y, double h) const;

    double estimateError(const StateVector& y_7, const StateVector& y_8) const;

    std::pair<StateVector, double> adaptiveStep(
        const RightHandSide& f, double t, const StateVector& y, double h, double h_min, double h_max) const;
};