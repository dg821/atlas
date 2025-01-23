//
// Created by Douglas Garza on 9/15/24.
//

#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>
#include "../../math/UniversalConstants.h"
#include "stateConversions.h"

namespace kepler {
    std::pair<Eigen::Vector3d, Eigen::Vector3d> solveKepler(const Eigen::Vector3d& r0, const Eigen::Vector3d& v0, double dt, double mu=UniversalConstants::EarthParams::MU);

    std::pair<double, double> findC2C3(double psi);
}

