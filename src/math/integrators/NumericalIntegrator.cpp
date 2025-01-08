//
// Created by Douglas Garza on 8/29/24.
//

#include "NumericalIntegrator.h"
#include <stdexcept>

NumericalIntegrator::NumericalIntegrator(StepType type)
    : stepType(type), atol(1e-12), rtol(1e-12) {}

void NumericalIntegrator::setTolerances(double atol, double rtol) {
    if (stepType != StepType::VariableStep) {
        throw std::logic_error("Tolerances are only applicable for variable step integrators");
    }
    this->atol = atol;
    this->rtol = rtol;
}

void NumericalIntegrator::setStepSizeLimits(double h_min, double h_max) {
    if (stepType != StepType::VariableStep) {
        throw std::logic_error("Step size limits are only applicable for variable step integrators");
    }
    minStepSize = h_min;
    maxStepSize = h_max;
}

NumericalIntegrator::StepType NumericalIntegrator::getStepType() const {
    return stepType;
}