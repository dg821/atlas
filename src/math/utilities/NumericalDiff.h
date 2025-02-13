//
// Created by Douglas Garza on 2/7/25.
//

#pragma once

#include <cmath>
#include <functional>
#include <stdexcept>

class NumericalDiff {
private:
    // Constants for adaptive step size
    static constexpr double EPSILON = 1e-14;  // Machine epsilon approximation
    static constexpr double MIN_STEP = 1e-12; // Minimum step size
    static constexpr int MAX_REFINEMENTS = 10;// Maximum step size refinements

    static double estimateStepSize(const std::function<double(double)>& f, double x) ;

    // Fifth-order central difference coefficients
    static constexpr int STENCIL_SIZE = 5;
    static constexpr double COEFFS[STENCIL_SIZE] = {
        1.0 / 12.0,   // h = ±2
        -2.0 / 3.0,   // h = ±1
        0.0,          // h = 0
        2.0 / 3.0,    // h = ±1
        -1.0 / 12.0   // h = ±2
    };


public:
    // Main differentiation function with error estimate
    static double differentiate( const std::function<double(double)>& f, double x, double& error_estimate, double tolerance = 1e-6);

    static double differentiate(const std::function<double(double)>& f, double x, double tolerance = 1e-6);
};