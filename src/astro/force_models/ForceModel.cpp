//
// Created by Douglas Garza on 8/25/24.
//

#include "ForceModel.h"

Eigen::Vector3d TwoBody::computeAcceleration(double t, const SpaceVehicle& sv) const {
    double mu = UniversalConstants::EarthParams::MU;
    SpaceVehicle::State state = sv.get_state().value();
    Eigen::Vector3d r = state.r;
    Eigen::Vector3d acceleration = (-mu / std::pow(r.norm(), 3)) * r;           // Two-body gravity model
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

    // Calculate acceleration components of J2 gravity model
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

    Eigen::Vector3d omegaEarthVec = {0.0, 0.0, omegaEarth};             // earth rotation vector

    Eigen::Vector3d relVel = v - omegaEarthVec.cross(r);                // relative velocity (accounts for earth rotation)
    double relVelMag = relVel.norm();
    double relVel2 = relVelMag * relVelMag;
    Eigen::Vector3d relVelUnit = relVel / relVelMag;

    AtmosphereDensityTable atmosphere;
    double hAboveEllipsoid = r.norm() - rEq;
    double airDensity = atmosphere.getDensity(hAboveEllipsoid) * 1e9;             // kg / km^3 (for unit agreement)

    Eigen::Vector3d acceleration = (-1.0/2.0) * (coeffDrag * area / mass) * airDensity * relVel2 * relVelUnit;      // drag model

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

    Eigen::Vector3d rSun = MeeusVectors::getSunVec(t);                  // sun vector from Meeus
    Eigen::Vector3d rSv2Sun = rSun - sv.get_state().value().r;

    Eigen::Vector3d acceleration = -(pSRP * reflectivity * area / mass) * (rSv2Sun / rSv2Sun.norm());

    return acceleration;
}

std::string SimpleSRPForce::getName() const {
    return "Simple SRP";
}

Eigen::Vector3d LunisolarThirdBodyForce::computeAcceleration(double t, const SpaceVehicle& sv) const {

    Eigen::Vector3d rSun = MeeusVectors::getSunVec(t);          // sun vector from Meeus
    Eigen::Vector3d rMoon = MeeusVectors::getMoonVec(t);        // moon vector from Meeus

    double muSun = UniversalConstants::SunParams::MU;
    double muMoon = UniversalConstants::MoonParams::MU;

    Eigen::Vector3d r = sv.get_state().value().r;                   // sv state

    Eigen::Vector3d rSv2Sun = rSun - r;                     // sv to sun
    Eigen::Vector3d rSv2Moon = rMoon - r;                   // sv to moon

    double rSunMag = rSun.norm();
    double rMoonMag = rMoon.norm();
    double rSv2SunMag = rSv2Sun.norm();
    double rSv2MoonMag = rSv2Moon.norm();

    Eigen::Vector3d accFromSun = muSun * (rSv2Sun / (rSv2SunMag * rSv2SunMag * rSv2SunMag) - rSun / (rSunMag * rSunMag * rSunMag));
    Eigen::Vector3d accFromMoon = muMoon * (rSv2Moon / (rSv2MoonMag * rSv2MoonMag * rSv2MoonMag) - rMoon / (rMoonMag * rMoonMag * rMoonMag));

    Eigen::Vector3d acceleration = accFromSun + accFromMoon;            // lunisolar third body model

    return acceleration;
}

std::string LunisolarThirdBodyForce::getName() const {
    return "Lunisolar Third Body";
}



void NonSphericalGravity::initializeCoefficients() {
    // Initialize with zeros
    Cnm = std::vector<std::vector<double>>(degree_n + 1,
        std::vector<double>(order_m + 1, 0.0));
    Snm = std::vector<std::vector<double>>(degree_n + 1,
        std::vector<double>(order_m + 1, 0.0));

    // load EGM96 normalized coefficients
    // Note: In production code, load complete set of coefficients
    Cnm[2][0] = -0.484165371736E-03;
    Cnm[2][2] = 0.243914352398E-05;
    Snm[2][2] = -0.140016683654E-05;
}

auto NonSphericalGravity::computeLegendrePolynomials(double phi) const ->
    std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> {

    std::vector<std::vector<double>> Pnm(degree_n + 1,
        std::vector<double>(order_m + 1, 0.0));
    std::vector<std::vector<double>> dPnm(degree_n + 1,
        std::vector<double>(order_m + 1, 0.0));

    double sin_phi = std::sin(phi);
    double cos_phi = std::cos(phi);
    double tan_phi = sin_phi/cos_phi;

    // seed values (n=0, m=0)
    Pnm[0][0] = 1.0;
    dPnm[0][0] = 0.0;

    // Compute P(1,0) and P(1,1) as base cases
    if (degree_n >= 1) {
        Pnm[1][0] = sqrt(3.0) * sin_phi;
        dPnm[1][0] = sqrt(3.0) * cos_phi;

        if (order_m >= 1) {
            Pnm[1][1] = sqrt(3.0) * cos_phi;
            dPnm[1][1] = -sqrt(3.0) * sin_phi;
        }
    }

    // compute polynomials using recursion relations
    for (int n = 2; n <= degree_n; n++) {
        // Compute P(n,0)
        double k = sqrt((2*n-1)/(2*n));
        Pnm[n][0] = k * sin_phi * Pnm[n-1][0];
        dPnm[n][0] = n * sin_phi * Pnm[n][0] - sqrt(n*(n+1)) * cos_phi * Pnm[n-1][0];

        // Compute P(n,m) and dP(n,m)
        for (int m = 1; m <= n && m <= order_m; m++) {
            if (n == m) {
                Pnm[n][m] = cos_phi * sqrt((2*n+1)/(2*n)) * Pnm[n-1][n-1];
                dPnm[n][m] = -n * tan_phi * Pnm[n][m];
            }
            else {
                double k1 = sqrt((2*n+1)/((2*n-1)*(n+m)*(n-m)));
                double k2 = sqrt((2*n+1)*(n+m-1)*(n-m-1)/((2*n-3)*(n+m)*(n-m)));
                Pnm[n][m] = k1 * sin_phi * Pnm[n-1][m] - k2 * Pnm[n-2][m];
                dPnm[n][m] = n * sin_phi * Pnm[n][m] - sqrt((n+m)*(n-m+1)) * cos_phi * Pnm[n][m-1];
            }
        }
    }

    return {Pnm, dPnm};
}

NonSphericalGravity::NonSphericalGravity()
    : NonSphericalGravity(DEFAULT_DEGREE, DEFAULT_ORDER) {  // delegate to the other constructor
}

NonSphericalGravity::NonSphericalGravity(int n, int m)
    : degree_n(std::min(n, MAX_DEGREE)), order_m(std::min(m, MAX_ORDER)) {
    if (n < 2 || m < 0 || m > n) {
        throw std::invalid_argument("Invalid gravity model orders");
    }
    initializeCoefficients();
}

Eigen::Vector3d NonSphericalGravity::computeAcceleration(double t, const SpaceVehicle& sv) const {
    double rEq = UniversalConstants::EarthParams::RADIUS_EQ;
    double mu = UniversalConstants::EarthParams::MU;

    // get position in ECI frame
    Eigen::Vector3d r_eci = sv.get_state().value().r;
    Eigen::Vector3d v_eci = sv.get_state().value().v;

    // transform from ECI to ECEF
    auto [r_ecef, v_ecef] = geodeticModel::eci2ecef(r_eci, v_eci, t);

    double x = r_ecef[0];
    double y = r_ecef[1];
    double z = r_ecef[2];

    double r = r_ecef.norm();
    double phi = std::asin(z/r);      // latitude
    double lambda = std::atan2(y, x);  // longitude

    // compute Legendre polynomials and derivatives
    auto [Pnm, dPnm] = computeLegendrePolynomials(phi);

    // initialize acceleration components (spherical coordinates)
    double ar = 0.0;    // radial
    double aphi = 0.0;  // latitudinal
    double alam = 0.0;  // longitudinal

    for (int n = 2; n <= degree_n; n++) {  // Start from n=2 as n=0,1 terms are not used
        double r_ratio = rEq/r;
        double r_ratio_n = std::pow(r_ratio, n);

        for (int m = 0; m <= n && m <= order_m; m++) {
            double cos_m_lambda = std::cos(m * lambda);
            double sin_m_lambda = std::sin(m * lambda);

            // radial component
            ar += (n+1) * r_ratio_n * Pnm[n][m] *
                 (Cnm[n][m] * cos_m_lambda + Snm[n][m] * sin_m_lambda);

            // latitudinal component
            aphi += r_ratio_n * dPnm[n][m] *
                   (Cnm[n][m] * cos_m_lambda + Snm[n][m] * sin_m_lambda);

            // longitudinal component
            if (m > 0) {
                alam += m * r_ratio_n * Pnm[n][m] *
                       (Snm[n][m] * cos_m_lambda - Cnm[n][m] * sin_m_lambda);
            }
        }
    }

    // Scale components
    ar *= -mu/(r*r);
    aphi *= -mu/(r*r);
    alam *= -mu/(r*r);

    // convert from spherical to ECEF
    double cos_phi = std::cos(phi);
    double sin_phi = std::sin(phi);
    double cos_lambda = std::cos(lambda);
    double sin_lambda = std::sin(lambda);

    // compute acceleration in ECEF
    Eigen::Vector3d a_ecef;
    a_ecef[0] = (ar * cos_phi * cos_lambda -
                     aphi * sin_phi * cos_lambda -
                     alam * sin_lambda);
    a_ecef[1] = (ar * cos_phi * sin_lambda -
                     aphi * sin_phi * sin_lambda +
                     alam * cos_lambda);
    a_ecef[2] = ar * sin_phi + aphi * cos_phi;

    // transform acceleration back to ECI frame
    auto [r_eci_out, v_eci_out, a_eci_out] = geodeticModel::ecef2eci(r_ecef, v_ecef, a_ecef, t);
    return a_eci_out;
}

void NonSphericalGravity::setOrders(int n, int m) {
    if (n < 2 || m < 0 || m > n || n > MAX_ORDER || m > MAX_ORDER) {
        throw std::invalid_argument("Invalid gravity model orders");
    }
    degree_n = n;
    order_m = m;
    initializeCoefficients();
}