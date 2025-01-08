//
// Created by Douglas Garza on 9/23/24.
//

#include "CoordinateFrames.h"

#include "astro/basic_astrodynamics/stateConversions.h"


Eigen::Matrix3d coordinateFrames::eci2rsw(Eigen::Vector3d& r_eci, Eigen::Vector3d& v_eci) {

    if (r_eci.cols() > 1) {
        r_eci.transpose();
    }
    if (v_eci.cols() > 1) {
        v_eci.transpose();
    }

    Eigen::Vector3d rHat = r_eci / r_eci.norm();
    Eigen::Vector3d wHat = r_eci.cross(v_eci) / ((r_eci.cross(v_eci)).norm());
    Eigen::Vector3d sHat = wHat.cross(rHat);

    Eigen::Matrix3d dcm_rsw2eci = Eigen::Matrix3d::Zero();
    dcm_rsw2eci << rHat, sHat, wHat;

    Eigen::Matrix3d dcm_eci2rsw = dcm_rsw2eci.transpose();

    return dcm_eci2rsw;
}

Eigen::Matrix3d coordinateFrames::eci2pqw(Eigen::Vector3d& r_eci, Eigen::Vector3d& v_eci) {

    if (r_eci.cols() > 1) {
        r_eci.transpose();
    }
    if (v_eci.cols() > 1) {
        v_eci.transpose();
    }

    double mu = UniversalConstants::EarthParams::MU;

    Eigen::Vector3d eccVec = ((v_eci.squaredNorm() - mu / r_eci.norm()) * r_eci - (r_eci.dot(v_eci)) * v_eci) / mu;

    Eigen::Vector3d pHat = eccVec / eccVec.norm();
    Eigen::Vector3d wHat = r_eci.cross(v_eci) / ((r_eci.cross(v_eci)).norm());
    Eigen::Vector3d qHat = wHat.cross(pHat);

    Eigen::Matrix3d dcm_pqw2eci = Eigen::Matrix3d::Zero();
    dcm_pqw2eci << pHat, qHat, wHat;

    Eigen::Matrix3d dcm_eci2pqw = dcm_pqw2eci.transpose();

    return dcm_eci2pqw;
}

Eigen::Matrix3d coordinateFrames::eci2ntw(Eigen::Vector3d& r_eci, Eigen::Vector3d& v_eci) {
    if (r_eci.cols() > 1) {
        r_eci.transpose();
    }
    if (v_eci.cols() > 1) {
        v_eci.transpose();
    }

    Eigen::Vector3d tHat = v_eci / v_eci.norm();
    Eigen::Vector3d wHat = r_eci.cross(v_eci) / ((r_eci.cross(v_eci)).norm());
    Eigen::Vector3d nHat = tHat.cross(wHat);

    Eigen::Matrix3d dcm_ntw2eci = Eigen::Matrix3d::Zero();
    dcm_ntw2eci << tHat, wHat, nHat;

    Eigen::Matrix3d dcm_eci2ntw = dcm_ntw2eci.transpose();

    return dcm_eci2ntw;
}

Eigen::Matrix3d eci2eqw(Eigen::Vector3d& r_eci, Eigen::Vector3d& v_eci, double mu = UniversalConstants::EarthParams::MU) {
    stateConversion::KeplerianElements kep = stateConversion::Cart2Kep(r_eci, v_eci, mu);
    Eigen::Matrix3d dcm_eci2eqw = RotationMatrices::rotZ(kep.node) * RotationMatrices::rotX(kep.inc) * RotationMatrices::rotZ(kep.node);

    return dcm_eci2eqw;
}

Eigen::Matrix3d coordinateFrames::ecf2sez(Eigen::Vector3d rSiteVector) {
    Eigen::Vector3d kHat = {0, 0, 1};

    Eigen::Vector3d zHat = rSiteVector / rSiteVector.norm();
    Eigen::Vector3d eHat = kHat.cross(zHat) / ((kHat.cross(zHat)).norm());
    Eigen::Vector3d sHat = eHat.cross(zHat);

    Eigen::Matrix3d dcm_sez2eci = Eigen::Matrix3d::Zero();
    dcm_sez2eci << sHat, eHat, zHat;

    Eigen::Matrix3d dcm_eci2sez = dcm_sez2eci.transpose();

    return dcm_eci2sez;
}

Eigen::Matrix3d pqw2rsw(double truA) {
    Eigen::Matrix3d dcm_pqw2rsw = RotationMatrices::rotZ(truA);

    return dcm_pqw2rsw;
}

Eigen::Matrix3d sez2body(double yaw, double roll, double pitch) {
    Eigen::Matrix3d dcm_sez2body = RotationMatrices::rotZ(yaw) * RotationMatrices::rotX(roll) * RotationMatrices::rotY(pitch);

    return dcm_sez2body;
}

