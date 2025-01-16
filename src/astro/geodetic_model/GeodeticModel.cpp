//
// Created by Douglas Garza on 10/1/24.
//

#include "GeodeticModel.h"

#include <numeric>


double geodeticModel::getParamC(double geodeticLatitude, double R = UniversalConstants::EarthParams::RADIUS_EQ, double Ecc = UniversalConstants::EarthParams::ECC_EARTH) {
    double paramC = R / std::sqrt(1 - (Ecc * Ecc) * std::pow(std::sin(geodeticLatitude), 2));
    return paramC;
}

double geodeticModel::getParamS(double geodeticLatitude, double R = UniversalConstants::EarthParams::RADIUS_EQ, double Ecc = UniversalConstants::EarthParams::ECC_EARTH) {
    double paramS = (R * (1 - Ecc * Ecc)) / std::sqrt(1 - (Ecc * Ecc) * std::pow(std::sin(geodeticLatitude), 2));
    return paramS;
}

double geodeticModel::getHorizontalR(double paramC, double heightEllipsoid, double geodeticLatitude) {
    double horizontalR = (paramC + heightEllipsoid) * std::cos(geodeticLatitude);
    return horizontalR;
}

double geodeticModel::getVerticalR(double paramS, double heightEllipsoid, double geodeticLatitude) {
    double verticalR = (paramS + heightEllipsoid) * std::sin(geodeticLatitude);
    return verticalR;
}

double geodeticModel::getRSiteMagnitude(double horizontalR, double verticalR) {
    double rSiteMagnitude = std::sqrt(horizontalR * horizontalR + verticalR * verticalR);
    return rSiteMagnitude;
}

double geodeticModel::geocentric2geodetic(double geocentricLatitude, double Ecc = UniversalConstants::EarthParams::ECC_EARTH) {
    double geodeticLatitude = atan(tan(geocentricLatitude) / (1 - Ecc * Ecc));
    return geodeticLatitude;
}

double geodeticModel::geodetic2geocentric(double geodeticLatitude, double Ecc = UniversalConstants::EarthParams::ECC_EARTH) {
    double geocentricLatitude = atan((1 - Ecc * Ecc) * tan(geodeticLatitude));
    return geocentricLatitude;
}

Eigen::Vector3d geodeticModel::lla2ecef(double geodeticLatitude, double longitude, double heightEllipsoid, double R = UniversalConstants::EarthParams::RADIUS_EQ, double Ecc = UniversalConstants::EarthParams::ECC_EARTH) {
    Eigen::Vector3d rECEF = Eigen::Vector3d::Zero();

    double paramC = getParamC(geodeticLatitude, R, Ecc);
    double paramS = getParamS(geodeticLatitude, R, Ecc);

    rECEF(0) = (paramC + heightEllipsoid) * std::cos(geodeticLatitude) * std::cos(longitude);
    rECEF(1) = (paramC + heightEllipsoid) * std::cos(geodeticLatitude) * std::sin(longitude);
    rECEF(2) = (paramS + heightEllipsoid) * std::sin(geodeticLatitude);

    return rECEF;
}

std::tuple<double, double, double> geodeticModel::ecef2lla(Eigen::Vector3d rECEF, double R_Eq = UniversalConstants::EarthParams::RADIUS_EQ, double R_Pol = UniversalConstants::EarthParams::RADIUS_POL, double Ecc = UniversalConstants::EarthParams::ECC_EARTH) {
    double x = rECEF(0);
    double y = rECEF(1);
    double z = rECEF(2);
    double a = R_Eq;
    double b = std::sqrt(R_Pol * (1 - Ecc * Ecc)) * MathFunctions::sign(z);
    double rVert = std::sqrt(x*x + y*y);

    double E = (b * z - (a*a - b*b)) / (a * rVert);

    double sinAlpha = y / rVert;
    double cosAlpha = x / rVert;
    double alpha = std::atan2(sinAlpha, cosAlpha);

    double longitude = alpha;

    double F = (b * z + (a*a - b*b)) / (a * rVert);
    double P = (4 * (E * F + 1)) / 3;
    double Q = 2 * (E*E - F*F);
    double D = P * P * P + Q * Q;

    double nu = 0.0;
    if (D > 0.0) {
        double nu = pow(std::sqrt(D) - Q, 1.0/3.0) - pow(std::sqrt(D) + Q, 1.0/3.0);
    } else if (D < 0.0) {
        double nu = 2 * std::sqrt(-P) * std::cos((1.0/3.0) * std::acos(Q / (P * std::sqrt(-P))));
    } else {
        std::exception('error: D is undefined.');
    }

    double G = 0.5 * (std::sqrt(E * E + nu) + E);
    double t = std::sqrt(G * G + ((F - nu*G) / (2 * G - E))) - G;

    double geodeticLatitude = atan((a * (1 - t * t)) / (2 * b * t));
    double heightEllipsoid = (rVert - a * t) * std::cos(geodeticLatitude) + (z - b) * std::sin(geodeticLatitude);

    return std::make_tuple(geodeticLatitude, longitude, heightEllipsoid);

}

double geodeticModel::getLocalSiderealTime(double timej2k, double longitude) {

    double jd = timeConversions::timej2k_to_jd(timej2k);
    double T_UT1 = (jd - 2451545.0) / 36525.0;
    long double GMST = 6731054841.0 + (876600.0 + 8640184812866.0) * T_UT1 + 0.093104 * T_UT1 * T_UT1 - 6.2e-6 * T_UT1 * T_UT1 * T_UT1;

    if (GMST > 86400.0) {
        while (GMST > 86400.0) {
            GMST -= 86400.0;
        }
    } else if (GMST < -86400.0) {
        while (GMST < -86400.0) {
            GMST += 86400.0;
        }
    }

    GMST = (GMST / 240) * UniversalConstants::D2R;

    if (GMST < 0) {
        GMST += (2 * std::numbers::pi);
    }

    double localSiderealTime = GMST + longitude;

    return localSiderealTime;
}


Eigen::Matrix3d geodeticModel::getDCM_itrf2tirs(double timej2k) {
    std::tuple tt_ut1_tuple = timeConversions::getAuxiliaryTerrestrialTime(timej2k);
    double t_tt = std::get<0>(tt_ut1_tuple);

    TimeAndPoleData tpdata;
    timeConversions::GregorianDate dateToday = tpdata.getCurrentUTCTime();

    // double timej2kToday = date_to_timej2k(dateToday);
    double x_p;
    double y_p;
    // if (timej2kToday > timej2k) {
    //     // WCAFTL: get x_p and y_p from database
    //     TimeAndPoleData tpdata;
    //     std::string poleDataURL = tpdata.fetchData(tpdata.getDeltaUt1andPoleUrl());
    //     std::vector<TimeAndPoleData::PoleData> poleData = tpdata.parseIERSPoleData(poleDataURL);
    //
    // } else {

    double mjd = timeConversions::jd_to_mjd(timeConversions::timej2k_to_jd(timej2k));

    double A = 2 * M_PI * (mjd-60593)/365.25;
    double C = 2 * M_PI * (mjd-60593)/435;

    x_p =  0.1234 + 0.0626 * std::cos(A) - 0.0587 * std::sin(A) + 0.0429 * std::cos(C) + 0.0955 * std::sin(C);
    y_p =  0.3606 - 0.0618 * std::cos(A) - 0.0557 * std::sin(A) + 0.0955 * std::cos(C) - 0.0429 * std::sin(C);
    // }

    double a_c = UniversalConstants::EarthParams::CHANDLER_WOBBLE;
    double a_a = UniversalConstants::EarthParams::ANNUAL_WOBBLE;

    double s_prime = -0.0015 * (a_c * a_c / 1.2 + a_a * a_a) * t_tt;

    x_p = x_p * UniversalConstants::ARCSECONDS_2_DEGREES * UniversalConstants::D2R;
    y_p = y_p * UniversalConstants::ARCSECONDS_2_DEGREES * UniversalConstants::D2R;
    s_prime = s_prime * UniversalConstants::ARCSECONDS_2_DEGREES * UniversalConstants::D2R;

    Eigen::Matrix3d dcm_itrf2tirs = RotationMatrices::rotZ(-s_prime) * RotationMatrices::rotY(x_p) * RotationMatrices::rotX(y_p);

    return dcm_itrf2tirs;
}


Eigen::Matrix3d geodeticModel::getDCM_tirs2cirs(double timej2k) {
    std::tuple tt_ut1_tuple = timeConversions::getAuxiliaryTerrestrialTime(timej2k);
    timeConversions::GregorianDate ut1_date = std::get<1>(tt_ut1_tuple);

    double jd_ut1 = timeConversions::date_to_jd(ut1_date.year, ut1_date.month, ut1_date.day, ut1_date.hour, ut1_date.minute, ut1_date.second);
    double earthRotationAngle = 2 * M_PI * (0.7790572732640 + 1.00273781191135448 * (jd_ut1 - 2451545.0));

    Eigen::Matrix3d dcm_tirs2cirs = RotationMatrices::rotZ(-earthRotationAngle);

    return dcm_tirs2cirs;
}

Eigen::Matrix3d geodeticModel::getDCM_cirs2gcrf(double timej2k) {
    // WCAFTL: include dX and dY from data

    double t_tt = timeConversions::getJulianCenturiesOfTT(timej2k);

    PrecessionNutation_XData PN_X;
    PrecessionNutation_YData PN_Y;
    PrecessionNutation_sData PN_s;

    constexpr std::array<double, 6> xPolyCoeffs = PN_X.polynomialCoeffs;
    constexpr std::array<double, 6> yPolyCoeffs = PN_Y.polynomialCoeffs;
    constexpr std::array<double, 6> sPolyCoeffs = PN_s.polynomialCoeffs;

    std::vector<auto> dataObject = {PN_X, PN_Y, PN_s};

    std::vector polyCoeffs = {xPolyCoeffs, yPolyCoeffs, sPolyCoeffs};
    std::array<double, 3> totalPolyAndNonPoly{};

    for (int i = 0; i < polyCoeffs.size(); i++) {
        double polyTotal = 0.0;
        double nonPolyTotal = 0.0;
        for (int j = 0; j < 6; j++) {
            polyTotal += polyCoeffs[i][j] * std::pow(t_tt, j);
        }

        std::array<double, 5> eachSummationTotal{};
        double summationRes = 0.0;
        std::array<double, 12> summationArray{};
        for (int j = 0; j < 5; j++) {
            eachSummationTotal[j] = 0.0;
            auto nonPolynomialCoeffs = std::get<j>(dataObject.nonPolynomialCoeffs);
            for (int k = 0; k < nonPolynomialCoeffs.size(); k++) {
                double sinTerm = nonPolynomialCoeffs[k][0];
                double cosTerm = nonPolynomialCoeffs[k][1];
                auto funcs = FundamentalArguments::getAllFundamentalArguments(t_tt);
                std::transform(nonPolynomialCoeffs[k][2], nonPolynomialCoeffs[k][14], funcs.begin(), summationArray.begin(), [](int a, int b) { return a * b; });
                summationRes = std::accumulate(summationArray.begin(), summationArray.end(), 0.0);
                eachSummationTotal[j] += (sinTerm * std::sin(summationRes) + cosTerm * std::cos(summationRes));
            }
            eachSummationTotal[j] = eachSummationTotal[j] * std::pow(t_tt, j);
            nonPolyTotal += eachSummationTotal[j];
        }

        totalPolyAndNonPoly[i] = polyTotal + nonPolyTotal;

    }

    double nodeMoon = FundamentalArguments::getNodeMoon(t_tt);
    double dXX = (-16617.0 - 1.6 * t_tt * t_tt + 0.7 * std::cos(nodeMoon)) * UniversalConstants::ARCSECONDS_2_DEGREES;
    double dYY = (-6819.2 - 141.9 * t_tt + 0.5 * std::sin(nodeMoon)) * UniversalConstants::ARCSECONDS_2_DEGREES;

    double XX = totalPolyAndNonPoly[0] + dXX;
    double YY = totalPolyAndNonPoly[1] + dYY;
    double ss = totalPolyAndNonPoly[2] - XX*YY/2;
    double aa = 0.5 + (1.0 / 8.0) * (XX * XX + YY * YY);

    Eigen::Matrix3d intermediateMatrix;
    intermediateMatrix << 1.0 - aa * XX * XX, -aa * XX * YY, XX,
                    -aa * XX * YY, 1. - aa * YY * YY, YY,
                    -XX, -YY, 1.0 - aa * (XX * XX + YY * YY);

    Eigen::Matrix3d dcm_cirs2gcrf = intermediateMatrix * RotationMatrices::rotZ(ss);

    return dcm_cirs2gcrf;
}

Eigen::Matrix3d geodeticModel::getDCM_itrf2gcrf(double timej2k) {
    Eigen::Matrix3d dcm_itrf2tirs = getDCM_itrf2tirs(timej2k);
    Eigen::Matrix3d dcm_tirs2cirs = getDCM_tirs2cirs(timej2k);
    Eigen::Matrix3d dcm_cirs2gcrf = getDCM_cirs2gcrf(timej2k);

    Eigen::Matrix3d dcm_itrf2gcrf = dcm_cirs2gcrf * dcm_tirs2cirs * dcm_itrf2tirs;

    return dcm_itrf2gcrf;
}


std::tuple<Eigen::Vector3d, Eigen::Vector3d> geodeticModel::ecef2eci(Eigen::Vector3d r_ecef, Eigen::Vector3d v_ecef, double timej2k) {
    // IAU-2000 CIO-Based (ITRF to GCRF)

    Eigen::Matrix3d dcm_itrf2gcrf = getDCM_itrf2gcrf(timej2k);
    Eigen::Vector3d r_eci = dcm_itrf2gcrf * r_ecef;

    Eigen::Matrix3d dcm_itrf2tirs = getDCM_itrf2tirs(timej2k);
    Eigen::Matrix3d dcm_tirs2cirs = getDCM_tirs2cirs(timej2k);
    Eigen::Matrix3d dcm_cirs2gcrf = getDCM_cirs2gcrf(timej2k);

    Eigen::Vector3d r_tirs = dcm_itrf2tirs * r_ecef;
    Eigen::Vector3d earthRotationVector = {0.0, 0.0, UniversalConstants::EarthParams::ROTATION_RATE};

    Eigen::Vector3d v_eci = dcm_cirs2gcrf * dcm_tirs2cirs * (dcm_itrf2tirs * v_ecef + earthRotationVector.cross(r_tirs));

    return std::make_tuple(r_eci, v_eci);
}

std::tuple<Eigen::Vector3d, Eigen::Vector3d> geodeticModel::eci2ecef(Eigen::Vector3d r_eci, Eigen::Vector3d v_eci, double timej2k) {
    Eigen::Matrix3d dcm_gcrf2itrf = getDCM_itrf2gcrf(timej2k).transpose();
    Eigen::Vector3d r_ecef = dcm_gcrf2itrf * r_eci;

    Eigen::Matrix3d dcm_itrf2tirs = getDCM_itrf2tirs(timej2k);
    Eigen::Matrix3d dcm_tirs2cirs = getDCM_tirs2cirs(timej2k);
    Eigen::Matrix3d dcm_cirs2gcrf = getDCM_cirs2gcrf(timej2k);

    Eigen::Vector3d r_tirs = dcm_itrf2tirs * r_ecef;
    Eigen::Vector3d earthRotationVector = {0.0, 0.0, UniversalConstants::EarthParams::ROTATION_RATE};

    Eigen::Vector3d v_ecef = dcm_itrf2tirs.transpose() * (dcm_tirs2cirs.transpose() * dcm_cirs2gcrf.transpose() * v_eci - earthRotationVector.cross(r_tirs));

    return std::make_tuple(r_ecef, v_ecef);
}



std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d> geodeticModel::ecef2eci(Eigen::Vector3d r_ecef, Eigen::Vector3d v_ecef, Eigen::Vector3d a_ecef, double timej2k) {
    // IAU-2000 CIO-Based (ITRF to GCRF)

    Eigen::Matrix3d dcm_itrf2gcrf = getDCM_itrf2gcrf(timej2k);
    Eigen::Vector3d r_eci = dcm_itrf2gcrf * r_ecef;

    Eigen::Matrix3d dcm_itrf2tirs = getDCM_itrf2tirs(timej2k);
    Eigen::Matrix3d dcm_tirs2cirs = getDCM_tirs2cirs(timej2k);
    Eigen::Matrix3d dcm_cirs2gcrf = getDCM_cirs2gcrf(timej2k);

    Eigen::Vector3d r_tirs = dcm_itrf2tirs * r_ecef;
    Eigen::Vector3d v_tirs = dcm_itrf2tirs * v_ecef;        // omega_polar is neglected here (assumed to be small)

    Eigen::Vector3d earthRotationVector = {0.0, 0.0, UniversalConstants::EarthParams::ROTATION_RATE};

    Eigen::Vector3d v_eci = dcm_cirs2gcrf * dcm_tirs2cirs * (dcm_itrf2tirs * v_ecef + earthRotationVector.cross(r_tirs));

    Eigen::Vector3d a_eci = dcm_cirs2gcrf * dcm_tirs2cirs * (dcm_itrf2tirs * a_ecef + 2 * earthRotationVector.cross(v_tirs) + earthRotationVector.cross(earthRotationVector.cross(r_tirs)));      // omegaDot term is neglected

    return std::make_tuple(r_eci, v_eci, a_eci);
}