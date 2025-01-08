//
// Created by Douglas Garza on 10/3/24.
//

#include "MathFunctions.h"


double sign(double x) {
    if (x > 0.0) {
        return 1.0;
    } else if (x < 0.0) {
        return -1.0;
    } else if (x == 0.0) {
        return 0.0;
    } else {
        std::exception('error: NaN.');
        return -2;
    }
}

inline double wrap2TwoPi(double angle) {
    double twoPi = 2.0 * M_PI;
    return angle - twoPi * floor( angle / twoPi );
}