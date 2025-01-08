//
// Created by Douglas Garza on 8/29/24.
//

#pragma once

#include <vector>
#include <functional>
#include <string>
#include <Eigen/Dense>

class NumericalIntegrator {
public:
    using Vector6d = Eigen::Matrix<double, 6, 1>;
    using StateVector = Vector6d;
    using RightHandSide = std::function<StateVector(double, const StateVector&)>;

    enum class StepType {
        FixedStep,
        VariableStep,
        MultiStep
    };

    explicit NumericalIntegrator(StepType type);
    virtual ~NumericalIntegrator() = default;

    // Pure virtual functions
    virtual StateVector step(const RightHandSide& f, double t, const StateVector& y, double h) const = 0;
    virtual std::vector<std::pair<double, StateVector>> integrate(const RightHandSide& f, double t_start, double t_end,
                                                            const StateVector& y_start, double h, double h_min, double h_max) const = 0;

    // Virtual functions with implementations
    virtual void setTolerances(double atol, double rtol);
    virtual void setStepSizeLimits(double h_min, double h_max);

    virtual int getOrder() const = 0;
    virtual std::string getName() const = 0;

    StepType getStepType() const;

protected:
    StepType stepType;
    double atol;
    double rtol;
    double minStepSize;
    double maxStepSize;
};