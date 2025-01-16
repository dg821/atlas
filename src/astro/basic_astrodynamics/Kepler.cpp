//
// Created by Douglas Garza on 9/15/24.
//

#include "Kepler.h"

namespace Kepler {
    std::pair<Eigen::Vector3d, Eigen::Vector3d> solveKepler(const Eigen::Vector3d& r0, const Eigen::Vector3d& v0, const double dt, const double mu) {
        constexpr double eps = 1e-8;            // convergence tolerance
        constexpr int maxIter = 50;             // maximum number of iterations

        const double r0_norm = r0.norm();
        const double v0_norm = v0.norm();

        // Calculate specific angular momentum
        const Eigen::Vector3d h = r0.cross(v0);
        double h_norm = h.norm();

        // Calculate eccentricity vector
        const Eigen::Vector3d eccVec = v0.cross(h) / mu - r0 / r0_norm;
        const double ecc = eccVec.norm();

        // Calculate semi-major axis
        const double sma = 1.0 / (2.0 / r0_norm - v0_norm * v0_norm / mu);

        // Calculate initial eccentric anomaly
        const double xi = r0.dot(v0) / std::sqrt(mu * sma);
        const double eta = r0_norm - sma;
        const double E0 = (xi * eta < 0) ? std::atan2(eta, xi) - M_PI : std::atan2(eta, xi);

        // Mean motion
        const double n = std::sqrt(mu / (sma * sma * sma));

        // Initial guess for E
        double E = E0 + n * dt;

        // Newton-Raphson iteration
        for (int iter = 0; iter < maxIter; ++iter) {
            const double f = E - ecc * std::sin(E) - n * dt;
            const double df = 1 - ecc * std::cos(E);
            const double dE = f / df;
            E -= dE;

            if (std::abs(dE) < eps) {
                break;
            }

            if (iter == maxIter - 1) {
                throw std::runtime_error("Newton-Raphson method did not converge");
            }
        }

        // Calculate new position and velocity in the orbital plane
        const double cosE = std::cos(E);
        const double sinE = std::sin(E);
        const double f = sma * (cosE - ecc);                                // f and g functions
        const double g = sma * std::sqrt(1 - ecc * ecc) * sinE;
        const double fdot = -sma * sma * n * sinE / r0_norm;
        const double gdot = sma * sma * n * (cosE - ecc) / r0_norm;

        // Calculate new position and velocity vectors
        Eigen::Vector3d r = f * r0 + g * v0;
        Eigen::Vector3d v = fdot * r0 + gdot * v0;

        return {r, v};
    }
}