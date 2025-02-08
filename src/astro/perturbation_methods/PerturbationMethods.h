//
// Created by Douglas Garza on 8/29/24.
//

#pragma once

#include "../../math/integrators/NumericalIntegrator.h"
#include "../../math/UniversalConstants.h"
#include "../SpaceVehicle/SpaceVehicle.h"
#include "../force_models/ForceModel.h"
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include "../basic_astrodynamics/kepler.h"
#include "../coordinate_frames/coordinateFrames.h"


class PerturbationMethod {
public:
    using Vector6d = Eigen::Matrix<double, 6, 1>;

    PerturbationMethod();
    virtual ~PerturbationMethod();
    virtual void updateState(SpaceVehicle& sv, double t);
    virtual std::string getName() const = 0;

    // Setters for integration method and force models
    void setIntegrationMethod(NumericalIntegrator* integrator);
    void setForceModels(const std::vector<ForceModel*>& models);

protected:
    NumericalIntegrator*integrator = nullptr;
    std::vector<ForceModel*> forceModels;
    TwoBody twoBodyForce;
    virtual Eigen::Vector3d computeTotalAcceleration(double t, Eigen::Vector3d& r, Eigen::Vector3d& v, SpaceVehicle& sv) const;
    virtual Eigen::Vector3d computePerturbAcceleration(double t, Eigen::Vector3d& r, Eigen::Vector3d& v, SpaceVehicle& sv) const;
    virtual Vector6d equationsOfMotion(double t, const Vector6d& stateVector, SpaceVehicle& sv) const = 0;
    double computeStepSize(double sma, double ecc, double mu) const;
    double computeMinStep(double sma, double ecc, double mu) const;
    double computeMaxStep(double sma, double ecc, double mu) const;
};

// Specific Perturbation Methods
class Cowell final : public PerturbationMethod {
public:
    std::string getName() const override;

    Vector6d equationsOfMotion(double t, const Vector6d& stateVector, SpaceVehicle& sv) const override;

    void updateState(SpaceVehicle& sv, double t) override;
};

// class Encke final : public PerturbationMethod {
//
// public:
//     void setEnckeStepSize(double dt) { m_dtEncke = dt; }
//     double getEnckeStepSize() const { return m_dtEncke.value_or(0.0); }
//
//     void setEnckeTolerance(double tol) { tolEncke = tol; }
//     double getEnckeTolerance() const { return tolEncke; }
//
//     std::string getName() const override;
//
//     Vector6d equationsOfMotion(double t, const Vector6d &stateVector, SpaceVehicle &sv, double f_val, double err, const Eigen::Vector3d &r_p, const
//                                Eigen::Vector3d &r) const ;
//
//     void updateState(SpaceVehicle& sv, double t) override;
//
// private:
//     std::optional<double> m_dtEncke;
//     double tolEncke = 0.01;
// };

// class VariationOfParameters final : public PerturbationMethod {
// public:
//     using Vector6d = Eigen::Matrix<double, 6, 1>;
//     std::string getName() const override;
//
//     Vector6d equationsOfMotionKeplerian(double t, const Vector6d& stateVector, SpaceVehicle& sv, double mu) const;
//     Vector6d equationsOfMotionEquinoctial(double t, const Vector6d& stateVector, SpaceVehicle& sv, double mu) const;
//     void updateState(SpaceVehicle& sv, double t) override;
// };