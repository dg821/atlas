//
// Created by Douglas Garza on 8/25/24.
//

#pragma once
#include <memory>
#include "../SpaceVehicle/SpaceVehicle.h"
#include "../force_models/ForceModel.h"
#include "../perturbation_methods/PerturbationMethods.h"
#include "../../math/integrators/NumericalIntegrator.h"
#include "../../math/integrators/Dopri87.h"

class OrbitPropagator {
public:

    OrbitPropagator();

    void setPerturbationMethod(std::unique_ptr<PerturbationMethod> method);
    void setIntegrationMethod(std::unique_ptr<NumericalIntegrator> method);
    void addForceModel(std::unique_ptr<ForceModel> model);

    void removeForceModel(const std::string& modelName);

    std::unique_ptr<PerturbationMethod> getPerturbationMethod();
    std::unique_ptr<NumericalIntegrator> getIntegrationMethod();
    std::vector<std::unique_ptr<ForceModel>> getForceModels();

    void propagate_to_time(SpaceVehicle& sv, double t);

    void propagate_to_perigee(SpaceVehicle& sv);
    void propagate_to_apogee(SpaceVehicle& sv);
    void propagate_to_angle(SpaceVehicle& sv, double theta);
    void propagate_to_sun(SpaceVehicle& sv);
    void propagate_to_moon(SpaceVehicle& sv);

private:
    void updatePerturbationMethodDependencies() const;

    std::unique_ptr<PerturbationMethod> perturbationMethod;
    std::unique_ptr<NumericalIntegrator> integrationMethod;
    std::vector<std::unique_ptr<ForceModel>> forceModels;
};
