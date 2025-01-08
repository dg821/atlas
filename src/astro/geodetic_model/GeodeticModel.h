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
#include "../../input_ouput/earth_time_and_pole_data/FundamentalArguments.h"
#include "../../input_ouput/earth_time_and_pole_data/PrecessionNutation_XData.h"
#include "../../input_ouput/earth_time_and_pole_data/PrecessionNutation_YData.h"
#include "../../input_ouput/earth_time_and_pole_data/PrecessionNutation_sData.h"

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
}