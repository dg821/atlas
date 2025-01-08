//
// Created by Douglas Garza on 1/7/25.
//

#pragma once

#include <vector>
#include <stdexcept>
#include <cmath>

struct AtmosphereLayer {
    double altitude_start;  // km
    double base_altitude;   // km
    double nominal_density; // kg/m^3
    double scale_height;    // km
};

class AtmosphereDensityTable {
private:
    std::vector<AtmosphereLayer> layers = {
        {0.0,   0,    1.225e0,    7.249},
        {25.0,  25,   3.899e-2,   6.349},
        {30.0,  30,   1.774e-2,   6.682},
        {40.0,  40,   3.972e-3,   7.554},
        {50.0,  50,   1.057e-3,   8.382},
        {60.0,  60,   3.206e-4,   7.714},
        {70.0,  70,   8.770e-5,   6.549},
        {80.0,  80,   1.905e-5,   5.799},
        {90.0,  90,   3.396e-6,   5.382},
        {100.0, 100,  5.297e-7,   5.877},
        {110.0, 110,  9.661e-8,   7.263},
        {120.0, 120,  2.438e-8,   9.473},
        {130.0, 130,  8.484e-9,   12.636},
        {140.0, 140,  3.845e-9,   16.149},
        {150.0, 150,  2.070e-9,   22.523},
        {180.0, 180,  5.464e-10,  29.740},
        {200.0, 200,  2.789e-10,  37.105},
        {250.0, 250,  7.248e-11,  45.546},
        {300.0, 300,  2.418e-11,  53.628},
        {350.0, 350,  9.518e-12,  53.298},
        {400.0, 400,  3.725e-12,  58.515},
        {450.0, 450,  1.585e-12,  60.828},
        {500.0, 500,  6.967e-13,  63.822},
        {600.0, 600,  1.454e-13,  71.835},
        {700.0, 700,  3.614e-14,  88.667},
        {800.0, 800,  1.170e-14,  124.64},
        {900.0, 900,  5.245e-15,  181.05},
        {1000.0, 1000, 3.019e-15,  268.00}
    };

public:
    // Get atmosphere layer for height above ellipsoid
    AtmosphereLayer getLayer(double altitudeAboveEllipsoid) const {
        if (altitudeAboveEllipsoid < 0.0 || altitudeAboveEllipsoid > 1000.0) {
            throw std::out_of_range("Altitude must be between 0 and 1000 km");
        }

        // Find the appropriate layer
        for (size_t i = 0; i < layers.size() - 1; i++) {
            if (altitudeAboveEllipsoid >= layers[i].altitude_start &&
                altitudeAboveEllipsoid < layers[i + 1].altitude_start) {
                return layers[i];
            }
        }

        // If we're at the last layer
        return layers.back();
    }

    // Calculate density at a specific altitude using exponential interpolation
    double getDensity(double altitudeAboveEllipsoid) const {
        AtmosphereLayer layer = getLayer(altitudeAboveEllipsoid);
        double delta_h = altitudeAboveEllipsoid - layer.base_altitude;
        return layer.nominal_density *
               std::exp(-delta_h / layer.scale_height);
    }
};


