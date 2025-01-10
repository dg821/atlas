//
// Created by Douglas Garza on 8/25/24.
//

#pragma once
#include <Eigen/Dense>
#include "../../math/UniversalConstants.h"
#include "../../math/utilities/MathFunctions.h"

namespace stateConversion {
    constexpr double mu = UniversalConstants::EarthParams::MU;

    struct KeplerianElements {
        double sma{}, ecc{}, inc{}, node{}, argP{}, truA{}, specE{}, semiparam{};
        Eigen::Vector3d ecc_vec, h_vec, n_vec;
    };

    struct CartesianState {
        Eigen::Vector3d r, v;
    };

    struct EquinoctialElements {
        double sma, a_f, a_g, h_e, k_e, eccLon;
    };

    double meanA2EccA(double meanA, double ecc);

    double truA2EccA(double truA, double ecc);

    double EccA2truA(double EccA, double ecc);

    KeplerianElements Cart2Kep(const Eigen::Vector3d& r, const Eigen::Vector3d& v, double mu = UniversalConstants::EarthParams::MU);

    CartesianState Kep2Cart(double sma, double ecc, double inc, double node, double argP, double truA, double mu = UniversalConstants::EarthParams::MU);

    EquinoctialElements Kep2Equinoctial(double sma, double ecc, double inc, double node, double argP, double truA);

    KeplerianElements Equinoctial2Kep(double sma, double a_f, double a_g, double h_e, double k_e, double eccLon);

}