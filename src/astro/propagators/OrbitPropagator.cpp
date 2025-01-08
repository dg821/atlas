//
// Created by Douglas Garza on 8/25/24.
//

#include "OrbitPropagator.h"

OrbitPropagator::OrbitPropagator() {
    // Set default perturbation method (e.g., Cowell's method)
    setPerturbationMethod(std::make_unique<Encke>());

    // Set default integration method (e.g., RK4)
    setIntegrationMethod(std::make_unique<Dopri87>());

    // Add default force model (e.g., Two-Body force)
    addForceModel(std::make_unique<TwoBody>());

    // Update dependencies to ensure everything is connected
    updatePerturbationMethodDependencies();
}

void OrbitPropagator::setPerturbationMethod(std::unique_ptr<PerturbationMethod> method) {
    perturbationMethod = std::move(method);
    updatePerturbationMethodDependencies();
}
void OrbitPropagator::setIntegrationMethod(std::unique_ptr<NumericalIntegrator> method) {
    integrationMethod = std::move(method);
    updatePerturbationMethodDependencies();
}

void OrbitPropagator::addForceModel(std::unique_ptr<ForceModel> model) {
    forceModels.push_back(std::move(model));
    updatePerturbationMethodDependencies();
}

void OrbitPropagator::removeForceModel(const std::string& modelName) {
    forceModels.erase(
        std::remove_if(forceModels.begin(), forceModels.end(),
            [&modelName](const auto& model) { return model->getName() == modelName; }),
        forceModels.end()
    );
    updatePerturbationMethodDependencies();
}

std::unique_ptr<PerturbationMethod> OrbitPropagator::getPerturbationMethod() {
    return std::move(perturbationMethod);
}

std::unique_ptr<NumericalIntegrator> OrbitPropagator::getIntegrationMethod() {
    return std::move(integrationMethod);
}

std::vector<std::unique_ptr<ForceModel>> OrbitPropagator::getForceModels() {
    return std::move(forceModels);
}

void OrbitPropagator::propagate_to_time(SpaceVehicle& sv, double t) {
    // Implement propagation logic using the force models
    // For each time step:
    //   1. Compute total acceleration from all force models
    //   2. Update spacecraft state
    if (perturbationMethod) {
        perturbationMethod->updateState(sv, t);
    } else {
        throw std::runtime_error("Perturbation method not set.");
    }
}

void OrbitPropagator::propagate_to_perigee(SpaceVehicle& sv) {

}

void OrbitPropagator::propagate_to_apogee(SpaceVehicle& sv) {

}

void OrbitPropagator::propagate_to_angle(SpaceVehicle& sv, double theta) {

}

void OrbitPropagator::propagate_to_sun(SpaceVehicle& sv) {

}

void OrbitPropagator::propagate_to_moon(SpaceVehicle& sv) {

}

void OrbitPropagator::updatePerturbationMethodDependencies() const {
    if (perturbationMethod) {
        perturbationMethod->setIntegrationMethod(integrationMethod.get());
        std::vector<ForceModel*> forceModelPtrs;
        for (const auto& model : forceModels) {
            forceModelPtrs.push_back(model.get());
        }
        perturbationMethod->setForceModels(forceModelPtrs);
    }
}
