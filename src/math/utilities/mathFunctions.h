//
// Created by Douglas Garza on 10/3/24.
//

#pragma once
#include <cmath>
#include <stdexcept>

namespace mathFunctions {
    double sign(double x);

    double wrap2TwoPi(double angle);

    double stableArcCos(double x);

    double factorial(int x);
}
