//
// Created by Douglas Garza on 8/25/24.
//

#include "stateConversions.h"

#include <cmath>
#include <Eigen/Dense>


double stateConversion::meanA2EccA(double meanA, double ecc) {

    meanA = mathFunctions::wrap2TwoPi(meanA);
    double tol = 1e-8;

    double EccA;
    if ((meanA < 0.0 && meanA > -M_PI) || meanA > M_PI) {
        double EccA = meanA - ecc;
    } else {
        double EccA = meanA + ecc;
    }

    double EccA_last = 10.0 * EccA;
    while (std::abs(EccA - EccA_last) < tol) {
        EccA = EccA_last + (meanA - EccA_last + ecc * std::sin(EccA_last)) / (1.0 - ecc * std::cos(EccA_last));
    }

    return EccA;
}

double stateConversion::truA2EccA(double truA, double ecc) {
    double EccA = std::atan2(std::sqrt(1 - ecc * ecc) * std::sin(truA), ecc + std::cos(truA));

    return EccA;
}

double stateConversion::EccA2truA(double EccA, double ecc) {

    double beta = ecc / (1 + std::sqrt(1 - ecc * ecc));
    double truA = EccA + 2.0 * std::atan2(beta * std::sin(EccA), 1.0 - beta * std::cos(EccA));

    return truA;
}


stateConversion::KeplerianElements Cart2Kep(const Eigen::Vector3d& r, const Eigen::Vector3d& v, const double mu) {
    stateConversion::KeplerianElements kep;

        kep.h_vec = r.cross(v);
        kep.ecc_vec = v.cross(kep.h_vec) / mu - r / r.norm();
        kep.ecc = kep.ecc_vec.norm();
        kep.n_vec = Eigen::Vector3d(0, 0, 1).cross(kep.h_vec);

        if (r.dot(v) > 0) {
            kep.truA = std::acos(r.dot(kep.ecc_vec) / (r.norm() * kep.ecc_vec.norm()));
        } else {
            kep.truA = -std::acos(r.dot(kep.ecc_vec) / (r.norm() * kep.ecc_vec.norm()));
        }

        kep.inc = std::acos(kep.h_vec(2) / kep.h_vec.norm());

        if (kep.n_vec.dot(Eigen::Vector3d(0, 1, 0)) > 0) {
            kep.node = std::acos(kep.n_vec.dot(Eigen::Vector3d(1, 0, 0)) / kep.n_vec.norm());
        } else {
            kep.node = -std::acos(kep.n_vec.dot(Eigen::Vector3d(1, 0, 0)) / kep.n_vec.norm());
        }

        if (kep.ecc_vec.dot(Eigen::Vector3d(0,0,1)) > 0) {
            kep.argP = std::acos(kep.n_vec.dot(kep.ecc_vec) / (kep.n_vec.norm() * kep.ecc_vec.norm()));
        } else {
            kep.argP = -std::acos(kep.n_vec.dot(kep.ecc_vec) / (kep.n_vec.norm() * kep.ecc_vec.norm()));
        }

        kep.specE = v.squaredNorm() / 2 - mu / r.norm();
        kep.sma = -mu / (2 * kep.specE);
        kep.semiparam = kep.h_vec.squaredNorm() / mu;

    return kep;
}


stateConversion::CartesianState Kep2Cart(const double sma, const double ecc, const double inc, const double node, const double argP, const double truA, const double mu) {
    stateConversion::CartesianState cart;

    const double r_norm = (sma*(1-std::pow(ecc,2)))/(1 + ecc*std::cos(truA));
    const double h_norm = std::sqrt(mu*sma*(1-std::pow(ecc,2)));
    const double semiparameter = sma*(1-std::pow(ecc,2));

    cart.r(0) = r_norm*(std::cos(node)*std::cos(argP+truA) - std::sin(node)*std::sin(argP+truA)*std::cos(inc));
    cart.r(1) = r_norm*(std::sin(node)*std::cos(argP+truA) + std::cos(node)*std::sin(argP+truA)*std::cos(inc));
    cart.r(2) = r_norm*(std::sin(inc)*std::sin(argP+truA));

    cart.v(0) = ((cart.r(0)*h_norm*ecc)/(r_norm*semiparameter))*std::sin(truA) - (h_norm/r_norm)*(std::cos(node)*std::sin(argP+truA) + std::sin(node)*std::cos(argP+truA)*std::cos(inc));
    cart.v(1) = ((cart.r(1)*h_norm*ecc)/(r_norm*semiparameter))*std::sin(truA) - (h_norm/r_norm)*(std::sin(node)*std::sin(argP+truA) - std::cos(node)*std::cos(argP+truA)*std::cos(inc));
    cart.v(2) = ((cart.r(2)*h_norm*ecc)/(r_norm*semiparameter))*std::sin(truA) + (h_norm/r_norm)*std::sin(inc)*std::cos(argP+truA);

    return cart;
}


stateConversion::EquinoctialElements Kep2Equinoctial(const double sma, const double ecc, const double inc, const double node, const double argP, const double truA) {
    stateConversion::EquinoctialElements equin;

    double eccAnom = stateConversion::truA2EccA(truA, ecc);

    equin.sma = sma;
    equin.a_f = ecc * std::cos(argP + node);
    equin.a_g = ecc * std::sin(argP + node);
    equin.h_e = std::tan(inc / 2) * std::cos(node);
    equin.k_e = std::tan(inc / 2) * std::sin(node);
    equin.eccLon = node + argP + eccAnom;

    return equin;
}

stateConversion::KeplerianElements Equinoctial2Kep(const double sma, const double a_f, const double a_g, const double h_e, const double k_e, const double eccLon) {
    stateConversion::KeplerianElements kep;

    kep.sma = sma;
    kep.ecc = std::sqrt(a_f * a_f + a_g * a_g);
    kep.inc = std::atan2(2.0 * std::sqrt(h_e * h_e + k_e * k_e), 1.0 - h_e * h_e - k_e * k_e);
    kep.node = std::atan2(k_e, h_e);
    kep.argP = std::atan2(a_g * h_e - a_f * k_e, a_f * h_e + a_g * k_e);

    double eccAnom = eccLon - kep.node - kep.argP;
    kep.truA = stateConversion::EccA2truA(eccAnom, kep.ecc);

    return kep;
}

