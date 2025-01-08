//
// Created by Douglas Garza on 8/24/24.
//

#pragma once
#include <Eigen/Dense>
#include "../basic_astrodynamics/stateConversions.h"
#include "../../math/UniversalConstants.h"
#include <cmath>

class SpaceVehicle {
  public:

  struct State {
    double timej2k;
    Eigen::Vector3d r;
    Eigen::Vector3d v;

    double sma;
    double ecc;
    double inc;
    double node;
    double argP;
    double truA;
  };

    explicit SpaceVehicle(double mass_total = 1000.0, double surface_area = 10.0, double mass_propellant = 400.0);

    void set_rv(double t, Eigen::Vector3d &r, Eigen::Vector3d &v);
    void set_orb(double t, double sma, double ecc, double inc, double node, double argP, double truA);
    void set_mass(double mass_total);
    void set_surface_area(double surface_area);
    void set_propellant_mass(double mass_propellant);

    [[nodiscard]] std::optional<State> get_state() const;
    [[nodiscard]] std::optional<double> get_time() const;
    [[nodiscard]] std::optional<double> get_mass() const;
    [[nodiscard]] std::optional<double> get_surface_area() const;
    [[nodiscard]] std::optional<double> get_propellant_mass() const;
    [[nodiscard]] double get_period(double) const;


  private:
    double mass_total_;
    double surface_area_;
    double mass_propellant_;
    std::optional<State> state_;
};
