//
// Created by Douglas Garza on 10/1/24.
//

#pragma once
#include <Eigen/Dense>
#include "../../math/utilities/RotationMatrices.h"
#include "../time_systems/TimeConversions.h"
#include <../../math/UniversalConstants.h>
#include <../time_systems/TimeConversions.h>
#include <cmath>
#include <numbers>
#include <tuple>
#include "../../math/utilities/MathFunctions.h"
#include "../../input_output/earth_time_and_pole_data/FundamentalArguments.h"
#include "../../input_output/earth_time_and_pole_data/PrecessionNutation_XData.h"
#include "../../input_output/earth_time_and_pole_data/PrecessionNutation_YData.h"
#include "../../input_output/earth_time_and_pole_data/PrecessionNutation_sData.h"

namespace geodeticModel {
    // earth model
    double geocentric2geodetic(double geocentricLatitude, double Ecc);
    double geodetic2geocentric(double geodeticLatitude, double Ecc);
    Eigen::Vector3d lla2ecef(double geodeticLatitude, double longitude, double heightEllipsoid, double R, double Ecc);
    std::tuple<double, double, double> ecef2lla(Eigen::Vector3d rECF, double R_Eq, double R_Pol, double Ecc);

    double getParamC(double geodeticLatitude, double R, double Ecc);
    double getParamS(double geodeticLatitude, double R, double Ecc);
    double getHorizontalR(double paramC, double heightEllipsoid, double geodeticLatitude);
    double getVerticalR(double paramS, double heightEllipsoid, double geodeticLatitude);
    double getRSiteMagnitude(double horizontalR, double verticalR);

    double getLocalSiderealTime(double timej2k, double longitude);

    Eigen::Matrix3d getDCM_itrf2tirs(double timej2k);
    Eigen::Matrix3d getDCM_tirs2cirs(double timej2k);
    Eigen::Matrix3d getDCM_cirs2gcrf(double timej2k);
    Eigen::Matrix3d getDCM_itrf2gcrf(double timej2k);

    std::tuple<Eigen::Vector3d, Eigen::Vector3d> geodeticModel::ecef2eci(Eigen::Vector3d r_ecef, Eigen::Vector3d v_ecef, double timej2k);
    std::tuple<Eigen::Vector3d, Eigen::Vector3d> eci2ecef(Eigen::Vector3d r_eci, Eigen::Vector3d v_eci, double timej2k);

    std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d> geodeticModel::ecef2eci(Eigen::Vector3d r_ecef, Eigen::Vector3d v_ecef, Eigen::Vector3d a_ecef, double timej2k);

}