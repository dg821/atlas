//
// Created by Douglas Garza on 9/23/24.
//

#pragma once
#include <Eigen/Dense>
#include <../../math/utilities/RotationMatrices.h>
#include <../geodetic_model/geodeticModel.h>
#include <../../math/UniversalConstants.h>

namespace coordinateFrames {

    // Satellite frame conversions
    Eigen::Matrix3d eci2rsw(Eigen::Vector3d& r_eci, Eigen::Vector3d& v_eci);
    Eigen::Matrix3d eci2pqw(Eigen::Vector3d& r_eci, Eigen::Vector3d& v_eci);
    Eigen::Matrix3d eci2ntw(Eigen::Vector3d& r_eci, Eigen::Vector3d& v_eci);
    Eigen::Matrix3d eci2eqw(Eigen::Vector3d& r_eci, Eigen::Vector3d& v_eci, double mu = UniversalConstants::EarthParams::MU);
    Eigen::Matrix3d pqw2rsw(double truA);

    // terrestrial frame conversions
    Eigen::Matrix3d ecf2sez(Eigen::Vector3d rSiteVector);

    Eigen::Matrix3d sez2body(double yaw, double roll, double pitch);
}