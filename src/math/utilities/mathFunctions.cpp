//
// Created by Douglas Garza on 10/3/24.
//

#include "mathFunctions.h"

namespace mathFunctions {
    double sign(double x) {
        if (x > 0.0) {
            return 1.0;
        } else if (x < 0.0) {
            return -1.0;
        } else if (x == 0.0) {
            return 0.0;
        } else {
            throw std::runtime_error("error: NaN.");
            return -2;
        }
    }

    double wrap2TwoPi(double angle) {
        double twoPi = 2.0 * M_PI;
        return angle - twoPi * floor( angle / twoPi );
    }

    double stableArcCos(double x) {
        if (x <= -1.0) return M_PI;
        if (x >= 1.0) return 0.0;
        return std::acos(x);
    }

    double factorial(int n) {
        if (n < 0) {
            throw std::invalid_argument("Factorial not defined for negative numbers");
        }
        if (n == 0) return 1.0;

        double result = 1.0;
        for (int i = 1; i <= n; i++) {
            result *= i;
        }
        return result;
    }
}