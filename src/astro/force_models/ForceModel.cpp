//
// Created by Douglas Garza on 8/25/24.
//

#include "ForceModel.h"

Eigen::Vector3d TwoBody::computeAcceleration(double t, const SpaceVehicle& sv) const {
    double mu = UniversalConstants::EarthParams::MU;
    SpaceVehicle::State state = sv.get_state().value();
    Eigen::Vector3d r = state.r;
    Eigen::Vector3d acceleration = (-mu / std::pow(r.norm(), 3)) * r;
    return acceleration;
}

std::string TwoBody::getName() const {
    return "TwoBody";
}

Eigen::Vector3d J2Force::computeAcceleration(double t, const SpaceVehicle& sv) const {

    double mu = UniversalConstants::EarthParams::MU;
    double rEq = UniversalConstants::EarthParams::RADIUS_EQ;
    double j2 = UniversalConstants::EarthParams::J2;

    SpaceVehicle::State state = sv.get_state().value();
    Eigen::Vector3d r = state.r;
    const double r2 = r.norm() * r.norm();
    const double r5 = r2 * r2 * r.norm();

    // Calculate z squared term
    const double z2 = r(2) * r(2);

    // Calculate common factor
    const double comFact = 1.5 * j2 * mu * (rEq * rEq) / r5;

    // Calculate acceleration components
    Eigen::Vector3d acceleration;
    acceleration(0) = comFact * r(0) * (5.0 * z2 / r2 - 1.0);
    acceleration(1) = comFact * r(1) * (5.0 * z2 / r2 - 1.0);
    acceleration(2) = comFact * r(2) * (5.0 * z2 / r2 - 3.0);

    return acceleration;
}

std::string J2Force::getName() const {
    return "J2";
}

Eigen::Vector3d DragExpForce::computeAcceleration(double t, const SpaceVehicle& sv) const {

    double omegaEarth = UniversalConstants::EarthParams::ROTATION_RATE;
    double rEq = UniversalConstants::EarthParams::RADIUS_EQ;

    SpaceVehicle::State state = sv.get_state().value();
    Eigen::Vector3d v = state.v;
    Eigen::Vector3d r = state.r;

    double mass = sv.get_mass().value();
    double area = sv.get_surface_area().value();

    double coeffDrag = 2.2;         // upper atmosphere flat plate model

    Eigen::Vector3d omegaEarthVec = {0.0, 0.0, omegaEarth};

    Eigen::Vector3d relVel = v - omegaEarthVec.cross(r);
    double relVelMag = relVel.norm();
    double relVel2 = relVelMag * relVelMag;
    Eigen::Vector3d relVelUnit = relVel / relVelMag;

    AtmosphereDensityTable atmosphere;
    double hAboveEllipsoid = r.norm() - rEq;
    double airDensity = atmosphere.getDensity(hAboveEllipsoid);

    Eigen::Vector3d acceleration = (-1.0/2.0) * (coeffDrag * area / mass) * airDensity * relVel2 * relVelUnit;

    return acceleration;
}

std::string DragExpForce::getName() const {
    return "DragExponential";
}



Eigen::Vector3d SimpleSRPForce::computeAcceleration(double t, const SpaceVehicle& sv) const {

    double pSRP = 4.57e-6;          // N/m^2

    Eigen::Vector3d acceleration;

    return acceleration;
}

std::string SimpleSRPForce::getName() const {
    return "Simple SRP";
}