//
// Created by Douglas Garza on 8/24/24.
//

// UniversalConstants.h

#pragma once
#include <numbers>

namespace UniversalConstants {

    // Mathematical constants
    constexpr double R2D = 180.0 / std::numbers::pi;
    constexpr double D2R = std::numbers::pi / 180.0;
    constexpr double ARCSECONDS_2_DEGREES = 1.0 / 3600.0;
    constexpr double ARCSECONDS_2_RAD = ARCSECONDS_2_DEGREES * D2R;
    constexpr double AU_2_KM = 149597870.7;  // 1 AU in kilometers

    // Physical constants
    constexpr double BIG_G = 6.67430e-20;  // Gravitational constant [km^3 kg^-1 s^-2]
    constexpr double SPEED_C = 299792.458;  // Speed of light [km/s]

    // Time constants
    constexpr int SECONDS_PER_MINUTE = 60;
    constexpr int SECONDS_PER_HOUR = 60;
    constexpr double UTC_2_TAI = 37.000000;
    constexpr double TAI_2_TT = 32.184;
    constexpr double UTC_2_GPS = UTC_2_TAI - 19.000000;

    // Earth-specific constants
    namespace EarthParams {
        constexpr double MU = 3.986004415e5;  // Earth's gravitational parameter [km^3/s^2]
        constexpr double RADIUS_EQ = 6378.1363;   // Earth's mean radius [km]
        constexpr double J2 = 1.08262668e-3;   // Earth J2
        constexpr double J3 = -2.5323e-6;   // Earth J3
        constexpr double J4 = -1.6204e-6;  // Earth J4
        constexpr double FLATTENING_FACTOR = 3.3528131e-3;  // Flattening factor []
        constexpr double ROTATION_RATE = 7.2921150e-5; // Earth's rotation rate [rad/s]
        constexpr double ECC_EARTH = 0.081819221456;    // Eccentricity []
        constexpr double RADIUS_POL = RADIUS_EQ * (1.0 - FLATTENING_FACTOR);
        constexpr double CHANDLER_WOBBLE = 0.26;        // Chandler Wobble [arcsec]
        constexpr double ANNUAL_WOBBLE = 0.12;          // Annual Wobble [arcsec]
    }

    // Sun-specific constants
    namespace SunParams {
        constexpr double MU = 1.32712428e11; // Sun's gravitational parameter [km^3/s^2]
        constexpr double RADIUS_EQ = 696000.0000;  // Sun's radius [km]
    }

    // Moon-specific constants
    namespace MoonParams {
        constexpr double MU = 4902.799;               // Moon's gravitational parameter [km^3/s^2]
        constexpr double RADIUS_EQ = 1738.0;            // Moon's radius [km]
    }
}
