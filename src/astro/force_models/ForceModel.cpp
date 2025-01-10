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

    double mass = sv.get_mass().value();                    // kg
    double area = sv.get_surface_area().value() * 1e-6;       // km (for unit agreement)

    double coeffDrag = 2.2;         // upper atmosphere flat plate model

    Eigen::Vector3d omegaEarthVec = {0.0, 0.0, omegaEarth};

    Eigen::Vector3d relVel = v - omegaEarthVec.cross(r);
    double relVelMag = relVel.norm();
    double relVel2 = relVelMag * relVelMag;
    Eigen::Vector3d relVelUnit = relVel / relVelMag;

    AtmosphereDensityTable atmosphere;
    double hAboveEllipsoid = r.norm() - rEq;
    double airDensity = atmosphere.getDensity(hAboveEllipsoid) * 1e9;             // kg / km^3 (for unit agreement)

    Eigen::Vector3d acceleration = (-1.0/2.0) * (coeffDrag * area / mass) * airDensity * relVel2 * relVelUnit;

    return acceleration;
}

std::string DragExpForce::getName() const {
    return "DragExponential";
}


Eigen::Vector3d SimpleSRPForce::computeAcceleration(double t, const SpaceVehicle& sv) const {

    double pSRP = 4.57e-6 * 1e3;          // N/m^2
    double reflectivity = 1.0;          // Approximated value. can vary between [0.0, 2.0]
    double area = sv.get_surface_area().value() * 1e-6;         // km (for unit agreement)
    double mass = sv.get_mass().value();

    Eigen::Vector3d rSun = MeeusVectors::getSunVec(t);
    Eigen::Vector3d rSv2Sun = rSun - sv.get_state().value().r;

    Eigen::Vector3d acceleration = -(pSRP * reflectivity * area / mass) * (rSv2Sun / rSv2Sun.norm());

    return acceleration;
}

std::string SimpleSRPForce::getName() const {
    return "Simple SRP";
}

Eigen::Vector3d LunisolarThirdBodyForce::computeAcceleration(double t, const SpaceVehicle& sv) const {

    Eigen::Vector3d rSun = MeeusVectors::getSunVec(t);
    Eigen::Vector3d rMoon = MeeusVectors::getMoonVec(t);

    double muSun = UniversalConstants::SunParams::MU;
    double muMoon = UniversalConstants::MoonParams::MU;

    Eigen::Vector3d r = sv.get_state().value().r;

    Eigen::Vector3d rSv2Sun = rSun - r;
    Eigen::Vector3d rSv2Moon = rMoon - r;

    double rSunMag = rSun.norm();
    double rMoonMag = rMoon.norm();
    double rSv2SunMag = rSv2Sun.norm();
    double rSv2MoonMag = rSv2Moon.norm();

    Eigen::Vector3d accFromSun = muSun * (rSv2Sun / (rSv2SunMag * rSv2SunMag * rSv2SunMag) - rSun / (rSunMag * rSunMag * rSunMag));
    Eigen::Vector3d accFromMoon = muMoon * (rSv2Moon / (rSv2MoonMag * rSv2MoonMag * rSv2MoonMag) - rMoon / (rMoonMag * rMoonMag * rMoonMag));

    Eigen::Vector3d acceleration = accFromSun + accFromMoon;

    return acceleration;
}

std::string LunisolarThirdBodyForce::getName() const {
    return "Lunisolar Third Body";
}