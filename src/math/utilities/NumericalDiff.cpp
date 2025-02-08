//
// Created by Douglas Garza on 2/7/25.
//

#include "NumericalDiff.h"

// Helper function to estimate optimal step size
double NumericalDiff::estimateStepSize(const std::function<double(double)>& f, double x) {
    double fx = std::abs(f(x));
    double h = (fx > EPSILON) ?
               std::sqrt(EPSILON) * std::max(std::abs(x), 1.0) :
               std::sqrt(EPSILON);
    return std::max(h, MIN_STEP);
}


// Main differentiation function with error estimate
double NumericalDiff::differentiate( const std::function<double(double)>& f, double x, double& error_estimate, double tolerance) {
    // Initial step size based on function and point characteristics
    double h = estimateStepSize(f, x);

    // Previous iteration results for convergence testing
    double prev_derivative = 0.0;
    double prev_error = INFINITY;

    // Refined calculation with step size adjustment
    for (int refinement = 0; refinement < MAX_REFINEMENTS; ++refinement) {
        // Calculate derivative using 5-point stencil
        double derivative = 0.0;
        for (int i = -2; i <= 2; ++i) {
            if (i != 0) {  // Skip center point
                derivative += COEFFS[i+2] * f(x + i*h);
            }
        }
        derivative /= h;

        // Calculate derivative with half step size for error estimation
        double half_h = h / 2.0;
        double half_derivative = 0.0;
        for (int i = -2; i <= 2; ++i) {
            if (i != 0) {  // Skip center point
                half_derivative += COEFFS[i+2] * f(x + i*half_h);
            }
        }
        half_derivative /= half_h;

        // Richardson extrapolation for error estimation
        error_estimate = std::abs(derivative - half_derivative) / 15.0;

        // Check for convergence
        if (error_estimate <= tolerance) {
            return derivative;
        }

        // Check if error is growing
        if (error_estimate > prev_error) {
            error_estimate = prev_error;
            return prev_derivative;
        }

        // Update for next iteration
        prev_derivative = derivative;
        prev_error = error_estimate;
        h /= 2.0;

        // Check if step size is too small
        if (h < MIN_STEP) {
            throw std::runtime_error("Failed to achieve desired accuracy: step size too small");
        }
    }

    throw std::runtime_error("Failed to converge within maximum refinements");
}

// Simplified interface without error estimate
double NumericalDiff::differentiate(const std::function<double(double)>& f, double x, double tolerance) {
    double error_estimate;
    return differentiate(f, x, error_estimate, tolerance);
}
