#pragma once

#include <cmath>
#include <array>
#include "../time_systems/TimeConversions.h"
#include "../../math/UniversalConstants.h"
#include <Eigen/Dense>

class MeeusVectors {
private:
    // calculate mean obliquity of the ecliptic (IAU 2006)
    static double getMeanObliquity(double T) {
        return (23.439291 - 0.0130042 * T - 1.64e-7 * T * T + 5.04e-7 * T * T * T) * UniversalConstants::D2R;
    }

    // transform from ecliptic to ECI coordinates
    static Eigen::Vector3d eclipticToECI(const Eigen::Vector3d& v_ecl, double T) {
        double eps = getMeanObliquity(T);
        double cos_eps = std::cos(eps);
        double sin_eps = std::sin(eps);

        return Eigen::Vector3d(
            v_ecl.x,
            v_ecl.y * cos_eps - v_ecl.z * sin_eps,
            v_ecl.y * sin_eps + v_ecl.z * cos_eps
        );
    }
public:
    // Meeus algorithm for Sun ephemeris simple analytic model. Accuracy 0.01 deg
    static Eigen::Vector3d getSunVec(double j2kseconds) {
        double T = timeConversions::secondsToT(j2kseconds);

        // Mean elements
        double L0 = 280.46646 + 36000.76983 * T + 0.0003032 * T * T;  // Mean longitude
        double M = 357.52911 + 35999.05029 * T - 0.0001537 * T * T;   // Mean anomaly
        double e = 0.016708634 - 0.000042037 * T - 0.0000001267 * T * T;  // Eccentricity

        // Reduce angles to 0-360 range
        L0 = std::fmod(L0, 360.0);
        M = std::fmod(M, 360.0);

        double M_rad = M * UniversalConstants::D2R;

        // Equation of center
        double C = (1.914602 - 0.004817 * T - 0.000014 * T * T) * std::sin(M_rad)
                + (0.019993 - 0.000101 * T) * std::sin(2.0 * M_rad)
                + 0.000289 * std::sin(3.0 * M_rad);

        // Sun's true longitude and true anomaly
        double L = L0 + C;
        double v = M + C;

        // Sun's radius vector in AU
        double R = (1.000001018 * (1.0 - e * e)) / (1.0 + e * std::cos(v * UniversalConstants::D2R));

        // Convert to rectangular coordinates
        double L_rad = L * UniversalConstants::D2R;
        Eigen::Vector3d sun_ecl = {
            R * std::cos(L_rad) * UniversalConstants::AU_2_KM,
            R * std::sin(L_rad) * UniversalConstants::AU_2_KM,
            0.0  // Assuming ecliptic coordinates
        };

        // Transform to ECI
        return eclipticToECI(sun_ecl, T);
    }

    // Meeus algorithm for Lunar ephemeris simple analytic model. Accuracy ~1 deg
    static Eigen::Vector3d getMoonVec(double j2kseconds) {
        double T = timeConversions::secondsToT(j2kseconds);

        // Lunar mean elements
        double Lp = 218.3164477 + 481267.88123421 * T  // Mean longitude
                  - 0.0015786 * T * T
                  + T * T * T / 538841.0
                  - T * T * T * T / 65194000.0;

        double D = 297.8501921 + 445267.1114034 * T    // Mean elongation
                  - 0.0018819 * T * T
                  + T * T * T / 545868.0
                  - T * T * T * T / 113065000.0;

        double M = 357.5291092 + 35999.0502909 * T     // Sun's mean anomaly
                  - 0.0001536 * T * T
                  + T * T * T / 24490000.0;

        double Mp = 134.9633964 + 477198.8675055 * T   // Moon's mean anomaly
                   + 0.0087414 * T * T
                   + T * T * T / 69699.0
                   - T * T * T * T / 14712000.0;

        double F = 93.2720950 + 483202.0175233 * T     // Argument of latitude
                  - 0.0036539 * T * T
                  - T * T * T / 3526000.0
                  + T * T * T * T / 863310000.0;

        // Reduce angles to 0-360 range
        Lp = std::fmod(Lp, 360.0);
        D = std::fmod(D, 360.0);
        M = std::fmod(M, 360.0);
        Mp = std::fmod(Mp, 360.0);
        F = std::fmod(F, 360.0);

        // Convert to radians
        double Lp_rad = Lp * UniversalConstants::D2R;
        double D_rad = D * UniversalConstants::D2R;
        double M_rad = M * UniversalConstants::D2R;
        double Mp_rad = Mp * UniversalConstants::D2R;
        double F_rad = F * UniversalConstants::D2R;

        // Calculate Moon's longitude
        double L = Lp + 6288.0160 * std::sin(Mp_rad)
                     + 1274.0097 * std::sin(2.0 * D_rad - Mp_rad)
                     + 658.7140 * std::sin(2.0 * D_rad)
                     + 214.2620 * std::sin(2.0 * Mp_rad)
                     - 186.0190 * std::sin(M_rad)
                     - 114.3320 * std::sin(2.0 * F_rad);

        // Calculate Moon's latitude
        double B = 5128.0 * std::sin(F_rad)
                + 280.0 * std::sin(Mp_rad + F_rad)
                + 277.0 * std::sin(Mp_rad - F_rad)
                + 176.0 * std::sin(2.0 * D_rad - F_rad)
                + 115.0 * std::sin(2.0 * D_rad + F_rad);

        // Calculate Moon's distance (in Earth radii)
        double R = 385000.56  // Moon's mean distance in km
                + 20905.355 * std::cos(Mp_rad)
                + 3699.111 * std::cos(2.0 * D_rad - Mp_rad)
                + 2955.968 * std::cos(2.0 * D_rad)
                + 569.925 * std::cos(2.0 * Mp_rad);

        // Convert to rectangular coordinates
        double L_rad = L * UniversalConstants::D2R;
        double B_rad = B * UniversalConstants::D2R;

        double cosB = std::cos(B_rad);
        Eigen::Vector3d moon_ecl = {
            R * cosB * std::cos(L_rad),
            R * cosB * std::sin(L_rad),
            R * std::sin(B_rad)
        };

        // Transform to ECI
        return eclipticToECI(moon_ecl, T);
    }
};