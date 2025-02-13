//
// Created by Douglas Garza on 8/25/24.
//

#include "ForceModel.h"

TwoBody::TwoBody() = default;
TwoBody::~TwoBody() = default;

Eigen::Vector3d TwoBody::computeAcceleration(double t, const Eigen::Vector3d& r, const Eigen::Vector3d& v, const SpaceVehicle& sv) const {

    double mu = UniversalConstants::EarthParams::MU;
    double rMag = r.norm();
    return (-mu / (rMag * rMag * rMag)) * r;           // Two-body gravity model

}

std::string TwoBody::getName() const {
    return "TwoBody";
}

DragExpForce::DragExpForce() = default;
DragExpForce::~DragExpForce() = default;

Eigen::Vector3d DragExpForce::computeAcceleration(double t, const Eigen::Vector3d& r, const Eigen::Vector3d& v, const SpaceVehicle& sv) const {

    double omegaEarth = UniversalConstants::EarthParams::ROTATION_RATE;
    double rEq = UniversalConstants::EarthParams::RADIUS_EQ;

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


SimpleSRPForce::SimpleSRPForce() = default;
SimpleSRPForce::~SimpleSRPForce() = default;

Eigen::Vector3d SimpleSRPForce::computeAcceleration(double t, const Eigen::Vector3d& r, const Eigen::Vector3d& v, const SpaceVehicle& sv) const {

    double pSRP = 4.57e-6 * 1e3;          // N/m^2
    double reflectivity = 1.0;          // Approximated value. can vary between [0.0, 2.0]
    double area = sv.get_surface_area().value() * 1e-6;         // km (for unit agreement)
    double mass = sv.get_mass().value();

    Eigen::Vector3d rSun = MeeusVectors::getSunVec(t);                  // sun vector from Meeus
    Eigen::Vector3d rSv2Sun = rSun - r;

    Eigen::Vector3d acceleration = -(pSRP * reflectivity * area / mass) * (rSv2Sun / rSv2Sun.norm());

    return acceleration;
}

std::string SimpleSRPForce::getName() const {
    return "Simple SRP";
}

LunisolarThirdBodyForce::LunisolarThirdBodyForce() = default;
LunisolarThirdBodyForce::~LunisolarThirdBodyForce() = default;

Eigen::Vector3d LunisolarThirdBodyForce::computeAcceleration(double t, const Eigen::Vector3d& r, const Eigen::Vector3d& v, const SpaceVehicle& sv) const {

    Eigen::Vector3d rSun = MeeusVectors::getSunVec(t);          // sun vector from Meeus
    Eigen::Vector3d rMoon = MeeusVectors::getMoonVec(t);        // moon vector from Meeus

    double muSun = UniversalConstants::SunParams::MU;
    double muMoon = UniversalConstants::MoonParams::MU;

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


NonSphericalGravity::~NonSphericalGravity() = default;

void NonSphericalGravity::initializeCoefficients() {
    // Initialize arrays with zeros
    Clm = std::vector<std::vector<double>>(degree_l + 1,
        std::vector<double>(order_m + 1, 0.0));
    Slm = std::vector<std::vector<double>>(degree_l + 1,
        std::vector<double>(order_m + 1, 0.0));

    if (!loadCoefficientsFromFile()) {
        throw std::runtime_error("Failed to load gravity coefficients from file: "
                               + coefficients_path);
    }
}

bool NonSphericalGravity::loadCoefficientsFromFile() {
    std::ifstream file(coefficients_path);
    if (!file.is_open()) {
        return false;
    }

    std::string line;
    // skip header
    std::getline(file, line);

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value;
        std::vector<double> row;

        // Parse CSV line
        while (std::getline(ss, value, ',')) {
            row.push_back(std::stod(value));
        }

        if (row.size() != 4) {  // l, m, C, S
            continue;
        }

        int l = static_cast<int>(row[0]);
        int m = static_cast<int>(row[1]);

        // Only store coefficients up to our maximum order
        if (l <= degree_l && m <= order_m) {
            Clm[l][m] = row[2];
            Slm[l][m] = row[3];
        }
    }

    return true;
}

double NonSphericalGravity::normalizeLegendreFunction(double Plm, int l, int m, int k) const {

    double Plm_normalized = std::sqrt((2 - k) * (2 * l + 1) *
        (mathFunctions::factorial(l - m) / mathFunctions::factorial(l + m))) * Plm;

    return Plm_normalized;
}

std::vector<std::vector<double>> NonSphericalGravity::computeLegendrePolynomials(double phi) const {
    // returns normalized legendre functions

    std::vector<std::vector<double>> Plm(degree_l + 1,
        std::vector<double>(degree_l + 1, 0.0));

    double sin_phi = std::sin(phi);
    double cos_phi = std::cos(phi);
    double sin_phi_squared = sin_phi * sin_phi;

    int l = degree_l;
    int m;
    int k;
    while (l >= 2) {

        double P_l_l_temp = mathFunctions::factorial(2 * l - 1) / (std::pow(2.0, l-1.0) *
            mathFunctions::factorial(l-1)) *
                std::pow(1.0 - sin_phi_squared, l/2.0);

        double P_l_lminus1_temp = sin_phi / std::sqrt(1.0 - sin_phi_squared) * P_l_l_temp;

        m = l;
        while (m >= 0) {
            if (m == l) {
                // if this loop just started
                if (m > 0) {
                    k = 0;
                } else {
                    k = 1;
                }
                Plm[l][l] = normalizeLegendreFunction(P_l_l_temp, l, l, k);
                Plm[l][l-1] = normalizeLegendreFunction(P_l_lminus1_temp, l, l-1, k);

                m -= 2;

            } else {

                double P_star_lm = 2.0 * (m + 1) * std::sqrt(1.0 / ((l + m + 1.0)*(l - m))) *
                    (sin_phi / std::sqrt(1.0 - sin_phi_squared)) * Plm[l][m + 1] -
                        std::sqrt((l + m + 2.0)*(l - m - 1.0)/((l + m + 1.0)*(l - m))) * Plm[l][m+2];

                if (m > 0) {
                    Plm[l][m] = P_star_lm;
                } else {
                    Plm[l][m] = (1.0 / std::sqrt(2.0)) * P_star_lm;
                }

                m--;
            }

        }

        l--;
    }

    return Plm;
}


NonSphericalGravity::NonSphericalGravity(int l, int m, const std::string& coeff_path)
    : degree_l(std::min(l, MAX_DEGREE)),
      order_m(m == -1 ? degree_l : std::min(m, MAX_ORDER)),
      coefficients_path(coeff_path) {

    if (l < 2 || (m != -1 && (m < 0 || m > l))) {
        throw std::invalid_argument("Invalid gravity model orders");
    }
    initializeCoefficients();
}

void NonSphericalGravity::setDegreeOrder(int l, int m) {
    if (l < 2 || m < 0 || m > l || l > MAX_DEGREE || m > MAX_ORDER) {
        throw std::invalid_argument("Invalid gravity model orders");
    }
    degree_l = l;
    order_m = m;
    initializeCoefficients();
}

std::string NonSphericalGravity::getName() const {
    return "NonSphericalGravity";
}

Eigen::Vector3d NonSphericalGravity::computeAcceleration(double t, const Eigen::Vector3d& r, const Eigen::Vector3d& v, const SpaceVehicle& sv) const {
    double rEq = UniversalConstants::EarthParams::RADIUS_EQ;
    double mu = UniversalConstants::EarthParams::MU;

    // transform from ECI to ECEF
    auto [r_ecef, v_ecef] = geodeticModel::eci2ecef(r, v, t);

    double x = r_ecef[0];
    double y = r_ecef[1];
    double z = r_ecef[2];

    double rMag = r_ecef.norm();
    double phi = std::asin(z/rMag);      // latitude
    double lambda = std::atan2(y, x);  // longitude

    double sinPhi = std::sin(phi);
    double cosPhi = std::cos(phi);

    // compute Legendre polynomials and derivatives
    auto Plm = computeLegendrePolynomials(phi);
    auto Plm_minus = computeLegendrePolynomials(phi - numDiffTol);
    auto Plm_plus = computeLegendrePolynomials(phi + numDiffTol);


    // initialize acceleration components (spherical coordinates)
    Eigen::Vector3d a_ecef;

    double dU_dr = 0.0;
    double dU_dphi = 0.0;
    double dU_dlambda = 0.0;

    for (int l = 2; l <= degree_l; l++) {
        double r_ratio = rEq/rMag;
        double r_ratio_l = std::pow(r_ratio, l);

        for (int m = 0; m <= l && m <= order_m; m++) {
            double cos_m_lambda = std::cos(m * lambda);
            double sin_m_lambda = std::sin(m * lambda);

            double Plm_deriv = (Plm_plus[l][m] - Plm_minus[l][m]) / (2.0 * numDiffTol);

            dU_dr += Plm[l][m] *
                (Clm[l][m] * cos_m_lambda + Slm[l][m] * sin_m_lambda);

            dU_dphi += Plm_deriv * (Clm[l][m] * cos_m_lambda + Slm[l][m] * sin_m_lambda);

            dU_dlambda += m * Plm[l][m] * (-Clm[l][m] * sin_m_lambda + Slm[l][m] * cos_m_lambda);

        }

        dU_dr *= (l + 1.0) * r_ratio_l;
        dU_dphi *= r_ratio_l;
        dU_dlambda *= r_ratio_l;
    }

    // Scale the derivatives and convert to acceleration
    dU_dr *= -mu/(rMag*rMag);
    dU_dphi *= mu/(rMag*rMag);
    dU_dlambda *= mu/(rMag*rMag);

    Eigen::Vector3d a_sph = {dU_dr, dU_dphi, dU_dlambda};
    Eigen::Matrix3d sph2car;

    // Precompute trigonometric values
    double cos_phi = std::cos(phi);
    double sin_phi = std::sin(phi);
    double cos_lambda = std::cos(lambda);
    double sin_lambda = std::sin(lambda);

    sph2car(0, 0) = cos_phi * cos_lambda;
    sph2car(0, 1) = -sin_phi * cos_lambda;
    sph2car(0, 2) = -sin_lambda;

    sph2car(1, 0) = cos_phi * sin_lambda;
    sph2car(1, 1) = -sin_phi * sin_lambda;
    sph2car(1, 2) = cos_lambda;

    sph2car(2, 0) = sin_phi;
    sph2car(2, 1) = cos_phi;
    sph2car(2, 2) = 0.0;

    a_ecef = sph2car * a_sph;

    Eigen::Matrix3d dcm_itrf2gcrf = geodeticModel::getDCM_itrf2gcrf(t);

    Eigen::Vector3d a_eci = dcm_itrf2gcrf * a_ecef;

    return a_eci;
}



//J2
Eigen::Vector3d J2::computeAcceleration(double t, const Eigen::Vector3d& r, const Eigen::Vector3d& v, const SpaceVehicle& sv) const {
    return nsg.computeAcceleration(t, r, v, sv);
}

std::string J2::getName() const {
    return "J2";
}

Eigen::Vector3d J2SimpleTest::computeAcceleration(double t, const Eigen::Vector3d& r, const Eigen::Vector3d& v, const SpaceVehicle& sv) const {
    double rEq = UniversalConstants::EarthParams::RADIUS_EQ;
    double mu = UniversalConstants::EarthParams::MU;
    double j2 = UniversalConstants::EarthParams::J2;

    // transform from ECI to ECEF
    auto [r_ecef, v_ecef] = geodeticModel::eci2ecef(r, v, t);

    double x = r_ecef[0];
    double y = r_ecef[1];
    double z = r_ecef[2];

    double rMag = r_ecef.norm();
    double phi = std::asin(z/rMag);      // latitude
    double lambda = std::atan2(y, x);  // longitude

    double sinPhi = std::sin(phi);
    double cosPhi = std::cos(phi);
    double tanPhi = sinPhi / cosPhi;

    // initialize acceleration components (spherical coordinates)
    Eigen::Vector3d a_ecef;

    a_ecef(0) = (-3 * j2 * mu * rEq * rEq * x / (2 * std::pow(rMag, 5))) * (1 - 5 * (z * z) / (rMag * rMag));
    a_ecef(1) = (-3 * j2 * mu * rEq * rEq * y / (2 * std::pow(rMag, 5))) * (1 - 5 * (z * z) / (rMag * rMag));
    a_ecef(2) = (-3 * j2 * mu * rEq * rEq * z / (2 * std::pow(rMag, 5))) * (3 - 5 * (z * z) / (rMag * rMag));

    Eigen::Matrix3d dcm_itrf2gcrf = geodeticModel::getDCM_itrf2gcrf(t);

    Eigen::Vector3d a_eci = dcm_itrf2gcrf * a_ecef;

    return a_eci;
}
std::string J2SimpleTest::getName() const {
    return "J2SimpleTest";
}