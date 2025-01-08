//
// Created by Douglas Garza on 8/24/24.
//

#include "SpaceVehicle.h"


SpaceVehicle::SpaceVehicle(const double mass_total, const double surface_area, const double mass_propellant)
    : mass_total_(mass_total), surface_area_(surface_area), mass_propellant_(mass_propellant) {}

void SpaceVehicle::set_rv(const double t, Eigen::Vector3d &r, Eigen::Vector3d &v) {

    // Function for conversion to orbital elements
    stateConversion::KeplerianElements kep = stateConversion::Cart2Kep(r, v);

    state_ = State{t, r, v, kep.sma, kep.ecc, kep.inc, kep.node, kep.argP, kep.truA};
}

void SpaceVehicle::set_orb(double t, double sma, double ecc, double inc, double node, double argP, double truA) {
    stateConversion::CartesianState cart = stateConversion::Kep2Cart(sma, ecc, inc, node, argP, truA);

    state_ = State{t, cart.r, cart.v, sma, ecc, inc, node, argP, truA};

}

void SpaceVehicle::set_mass(const double mass_total) {
    mass_total_ = mass_total;
}

void SpaceVehicle::set_surface_area(const double surface_area) {
    surface_area_ = surface_area;
}

void SpaceVehicle::set_propellant_mass(double mass_propellant) {
    mass_propellant_ = mass_propellant;
}

std::optional<SpaceVehicle::State> SpaceVehicle::get_state() const {
    return state_;
}

std::optional<double> SpaceVehicle::get_time() const {
    return state_->timej2k;
}

std::optional<double> SpaceVehicle::get_mass() const {
    return mass_total_;
}

std::optional<double> SpaceVehicle::get_surface_area() const {
    return surface_area_;
}

std::optional<double> SpaceVehicle::get_propellant_mass() const {
    return mass_propellant_;
}

double SpaceVehicle::get_period(double mu) const {
    double sma = state_->sma;
    double period = 2 * M_PI * std::sqrt(std::pow(sma, 3) / mu);

    return period;
}
