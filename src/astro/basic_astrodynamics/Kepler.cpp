//
// Created by Douglas Garza on 9/15/24.
//

#include "Kepler.h"

namespace Kepler {
    std::pair<Eigen::Vector3d, Eigen::Vector3d> solveKepler(const Eigen::Vector3d& r0, const Eigen::Vector3d& v0, const double dt, const double mu) {
        constexpr int maxIter = 50;             // maximum number of iterations

        const double r0_norm = r0.norm();
        const double v0_norm = v0.norm();

        // Calculate initial eccentric anomaly
        const double alpha = -(v0_norm * v0_norm) / mu + 2.0 / r0_norm;
        double psi;
        double c2;
        double c3;
        double rCalc;

        // xi initial guess
        if (std::abs(alpha) < 1e-6) {
            // parabola or hyperbola
            throw std::runtime_error("Currently no support for parabola/hyperbola.");
        }
        double xi = std::sqrt(mu) * dt * alpha;

        double xiLast = -1.0 / 0.0;
        int iter = 0;
        while (std::abs(xi - xiLast) >= 1e-6) {
            xiLast = xi;

            psi = xi * xi * alpha;
            auto c2c3_res = findC2C3(psi);
            c2 = c2c3_res.first;
            c3 = c2c3_res.second;

            rCalc = xi * xi * c2 + (r0.dot(v0)/std::sqrt(mu)) * xi * (1.0 - psi*c3) + r0_norm * (1.0 - psi * c2);
            xi = xi + (std::sqrt(mu) * dt - xi * xi * xi * c3 - (r0.dot(v0) / std::sqrt(mu)) * xi * xi * c2 - r0_norm * xi * (1.0 - psi * c3)) / rCalc;

            iter++;
            if (iter == maxIter) {
                throw std::runtime_error("Newton-Raphson method did not converge");
            }
        }

        double f = 1.0 - (xi * xi / r0_norm) * c2;
        double g = dt - (xi * xi * xi / std::sqrt(mu)) * c3;
        double gDot = 1.0 - (xi * xi / rCalc) * c2;
        double fDot = (std::sqrt(mu) / (rCalc * r0_norm)) * xi * (psi * c3 - 1.0);

        Eigen::Vector3d r = f * r0 + g * v0;
        Eigen::Vector3d v = fDot * r0 + gDot * v0;

        double fgEval = f * gDot - fDot * g;
        if (fgEval - 1.0 > 1.0e-5) {
            throw std::runtime_error("F and G function error.");
        }

        return {r, v};
    }

    std::pair<double, double> findC2C3(const double psi) {

        double c2;
        double c3;
        double sqrtPsi = std::sqrt(psi);
        double sqrtNegPsi = std::sqrt(-psi);
        if (psi > 1e-6) {
            c2 = (1.0 - std::cos(sqrtPsi)) / psi;
            c3 = (sqrtPsi - std::sin(sqrtPsi)) / std::sqrt(psi * psi * psi);
        } else {
            if (psi < -1e-6) {
                c2 = (1.0 - std::cosh(sqrtNegPsi)) / psi;
                c3 = (std::sinh(sqrtNegPsi) - sqrtNegPsi) / std::sqrt(sqrtNegPsi * sqrtNegPsi * sqrtNegPsi);
            } else {
                c2 = 0.5;
                c3 = 1.0 / 6.0;
            }
        }

        return {c2, c3};
    }
}
