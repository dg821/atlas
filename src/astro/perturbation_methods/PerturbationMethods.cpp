//
// Created by Douglas Garza on 8/29/24.
//

#include "PerturbationMethods.h"

PerturbationMethod::PerturbationMethod() {}

PerturbationMethod::~PerturbationMethod() = default;
void PerturbationMethod::updateState(SpaceVehicle& sv, double t) {

}

// Setters for integration method and force models
void PerturbationMethod::setIntegrationMethod(NumericalIntegrator* integrator) { this->integrator = integrator; }
void PerturbationMethod::setForceModels(const std::vector<ForceModel*>& models) { this->forceModels = models; }

// Step sizes
double PerturbationMethod::computeStepSize(const double sma, const double ecc, const double mu = UniversalConstants::EarthParams::MU) const {
    const double period = 2.0 * M_PI * std::sqrt(std::pow(sma, 3) / mu);
    const double fraction = std::pow(10, log10(sma) - 1);
    return (period / fraction) * std::pow(1.0 - ecc, 2);
}

double PerturbationMethod::computeMinStep(const double sma, const double ecc, const double mu = UniversalConstants::EarthParams::MU) const {
    const double period = 2.0 * M_PI * std::sqrt(std::pow(sma, 3) / mu);
    const double fraction = std::pow(10, log10(sma) + 2);
    return (period / fraction) * std::pow(1.0 - ecc, 3);
}

double PerturbationMethod::computeMaxStep(const double sma, const double ecc, const double mu = UniversalConstants::EarthParams::MU) const {
    const double period = 2.0 * M_PI * std::sqrt(std::pow(sma, 3) / mu);
    const double fraction = std::pow(10, floor(log10(sma) - 1));
    return (period / fraction) * std::pow(1.0 - ecc, 2);
}

Eigen::Vector3d PerturbationMethod::computeTotalAcceleration(double t, Eigen::Vector3d& r, Eigen::Vector3d& v, SpaceVehicle& sv) const {
    Eigen::Vector3d totalAcceleration = Eigen::Vector3d::Zero();
    for (auto forceModel : forceModels) {
        totalAcceleration = totalAcceleration + forceModel->computeAcceleration(t, r, v, sv);
    }
    return totalAcceleration;
}


Eigen::Vector3d PerturbationMethod::computePerturbAcceleration(double t, Eigen::Vector3d& r, Eigen::Vector3d& v, SpaceVehicle& sv) const {
    Eigen::Vector3d totalAcceleration = Eigen::Vector3d::Zero();
    for (const auto forceModel : forceModels) {
        if (forceModel->getName() == "TwoBody") {
            continue;
        }
        totalAcceleration = totalAcceleration + forceModel->computeAcceleration(t, r, v, sv);
    }
    return totalAcceleration;
}

std::string Cowell::getName() const {
    return "Cowell";
}

PerturbationMethod::Vector6d Cowell::equationsOfMotion(const double t, const Vector6d& stateVector, SpaceVehicle& sv) const {
    // Extract position and velocity from state vector
    Eigen::Vector3d r(stateVector(0), stateVector(1), stateVector(2));
    Eigen::Vector3d v(stateVector(3), stateVector(4), stateVector(5));

    // Compute total acceleration from force models
    Eigen::Vector3d pertAcceleration = computePerturbAcceleration(t, r, v, sv);
    Eigen::Vector3d totalAcceleration = computeTotalAcceleration(t, r, v, sv);

    // Create and return the derivative of the state (dydt)
    Vector6d dydt(6);
    dydt << v(0), v(1), v(2),  // dx/dt = velocity
            totalAcceleration(0), totalAcceleration(1), totalAcceleration(2);  // dv/dt = acceleration
    return dydt;
}

// Specific Perturbation Methods
void Cowell::updateState(SpaceVehicle& sv, const double t) {
    SpaceVehicle::State svState = sv.get_state().value();
    Eigen::Vector3d r = svState.r;
    Eigen::Vector3d v = svState.v;
    Vector6d stateVector(6);
    stateVector << r(0), r(1), r(2), v(0), v(1), v(2);

    const double sma = svState.sma;
    const double ecc = svState.ecc;

    const double dt = computeStepSize(sma, ecc);
    const double dt_max = computeMaxStep(sma, ecc);
    const double dt_min = computeMinStep(sma, ecc);

    const double tStart = sv.get_time().value();
    const double tEnd = t;

    auto solution = integrator->integrate(
    [this, &sv](const double t, const Vector6d& stateVector) -> Eigen::VectorXd {
        return this->equationsOfMotion(t, stateVector, sv);
        },
    tStart, tEnd, stateVector, dt, dt_min, dt_max
    );

    auto newState = solution.back();

    Eigen::Vector3d new_r = {newState.second(0), newState.second(1), newState.second(2)};
    Eigen::Vector3d new_v = {newState.second(3), newState.second(4), newState.second(5)};

    sv.set_rv(t, new_r, new_v);
}





// std::string Encke::getName() const {
//     return "Encke";
// }
//
// PerturbationMethod::Vector6d Encke::equationsOfMotion(const double t, const Vector6d& stateVector, SpaceVehicle& sv, const double f_val, const double err, const Eigen::Vector3d& r_p, const Eigen::Vector3d& r) const {
//     double mu = UniversalConstants::EarthParams::MU;
//
//     // Extract position and velocity from state vector
//     Eigen::Vector3d del_r = stateVector.head<3>();
//     Eigen::Vector3d del_v = stateVector.tail<3>();
//
//     // Compute total acceleration from force models
//     Eigen::Vector3d perturbedAcceleration = computePerturbAcceleration(t, sv);
//     Eigen::Vector3d rddot = perturbedAcceleration + (mu / r.norm()) * (f_val * err * r_p - del_r);
//
//     // Create and return the derivative of the state (dydt)
//     Vector6d dydt;
//     dydt << del_v, rddot;
//     return dydt;
// }
//
// void Encke::updateState(SpaceVehicle& sv, double t) {
//     using Vector6d = Eigen::Matrix<double, 6, 1>;
//
//     const SpaceVehicle::State& svState = sv.get_state().value();
//     Eigen::Vector3d r = svState.r;
//     Eigen::Vector3d v = svState.v;
//
//     const double sma = svState.sma;
//     const double ecc = svState.ecc;
//
//     const double dt = computeStepSize(sma, ecc);
//     const double dt_max = computeMaxStep(sma, ecc);
//     const double dt_min = computeMinStep(sma, ecc);
//
//     const double dtEncke = m_dtEncke.value_or(dt_max);
//
//     const double tStart = sv.get_time().value();
//     const double tEnd = t;
//
//     double tNow = tStart;
//     Vector6d del_stateVector = Vector6d::Zero();
//     Eigen::Vector3d r_p = r;
//     Eigen::Vector3d v_p = v;
//
//     while (tNow < tEnd) {
//         const auto [r_twoBody, v_twoBody] = kepler::solveKepler(r, v, tNow);
//
//         const Eigen::Vector3d del_r = del_stateVector.head<3>();
//         const double r_twoBody_mag = r_twoBody.norm();
//         const double err = r_twoBody.dot(del_r) / (r_twoBody_mag * r_twoBody_mag);
//         const double f_val = (1.0 / err) * (1.0 - (1.0 / std::pow(1.0 - 2.0 * err, 1.5)));
//
//         integrator->integrate(
//             [this, &sv, &r_p, &r_twoBody, f_val, err](double t, const Vector6d& del_stateVector) -> Vector6d {
//                 return this->equationsOfMotion(t, del_stateVector, sv, f_val, err, r_p, r_twoBody);
//             },
//             tNow, std::min(tNow + dtEncke, tEnd), del_stateVector, dt, dt_min, dt_max
//         );
//
//         const Eigen::Vector3d new_del_r = del_stateVector.head<3>();
//         const Eigen::Vector3d new_del_v = del_stateVector.tail<3>();
//
//         if (new_del_r.norm() / r_p.norm() > tolEncke) {
//             r = r_p;
//             v = v_p;
//             del_stateVector.setZero();
//         } else {
//             r_p = r_twoBody + new_del_r;
//             v_p = v_twoBody + new_del_v;
//             tNow = std::min(tNow + dtEncke, tEnd);
//         }
//     }
//
//     sv.set_rv(tEnd, r_p, v_p);
// }


// std::string VariationOfParameters::getName() const {
//     return "VariationOfParameters";
// }
//
// PerturbationMethod::Vector6d VariationOfParameters::equationsOfMotionKeplerian(const double t, const Vector6d& stateVector, SpaceVehicle& sv, double mu = UniversalConstants::EarthParams::MU) const {
//     double sma = stateVector(0);
//     double ecc = stateVector(1);
//     double inc = stateVector(2);
//     double node = stateVector(3);
//     double argP = stateVector(4);
//     double truA = stateVector(5);
//
//     Eigen::Vector3d r = sv.get_state().value().r;
//     Eigen::Vector3d v_eci = sv.get_state().value().v;
//     Eigen::Matrix3d dcm_eci2rsw = coordinateFrames::eci2rsw( r_eci,  v_eci);
//
//     // Compute total acceleration from force models
//     Eigen::Vector3d perturbedAcceleration_eci = computePerturbAcceleration(t,  sv);
//     Eigen::Vector3d perturbedAcceleration_rsw = dcm_eci2rsw * perturbedAcceleration_eci;
//
//     double F_R = perturbedAcceleration_rsw(0);
//     double F_S = perturbedAcceleration_rsw(1);
//     double F_W = perturbedAcceleration_rsw(2);
//
//     double rMag = r_eci.norm();
//     double n = std::sqrt(mu / (sma * sma * sma));
//     double p = sma * (1.0 - ecc * ecc);
//     double argLat = argP + truA;
//     Eigen::Vector3d hVec = r_eci.cross(v_eci);
//     double hMag = hVec.norm();
//
//     double da_dt = (2.0 / (n * std::sqrt(1.0 - ecc*ecc))) * (ecc * std::sin(truA) * F_R + (p / rMag) * F_S);
//     double de_dt = (std::sqrt(1.0 - ecc*ecc) / (n * sma)) * (std::sin(truA) * F_R + (std::cos(truA) + (ecc + std::cos(truA)) / (1.0 + ecc * std::cos(truA))) * F_S);
//     double di_dt = ((rMag * std::cos(argLat)) / (n * sma*sma * std::sqrt(1.0 - ecc*ecc))) * F_W;
//     double dnode_dt = ((rMag * std::sin(argLat)) / (n * sma*sma * std::sqrt(1.0 - ecc*ecc) * std::sin(inc))) * F_W;
//     double dargP_dt = (std::sqrt(1.0 - ecc*ecc) / (n * sma * ecc)) * (-std::cos(truA) * F_R + std::sin(truA) * (1.0 + rMag/p) * F_S)
//             - ((rMag * std::sin(argLat)) / (hMag * std::tan(inc))) * F_W;
//     double dtruA_dt = hMag / (rMag * rMag) + (1.0 / (ecc * hMag)) * (p * std::cos(truA)) * F_R - (p + rMag) * std::sin(truA) * F_S;
//
//     // Create and return the derivative of the state (dydt)
//     Vector6d dydt;
//     dydt << da_dt, de_dt, di_dt, dnode_dt, dargP_dt, dtruA_dt;
//     return dydt;
// }
//
//
// PerturbationMethod::Vector6d VariationOfParameters::equationsOfMotionEquinoctial(const double t, const Vector6d& stateVector, SpaceVehicle& sv, double mu = UniversalConstants::EarthParams::MU) const {
//     double sma = stateVector(0);
//     double a_f = stateVector(1);
//     double a_g = stateVector(2);
//     double h_e = stateVector(3);
//     double k_e = stateVector(4);
//     double eccLon = stateVector(5);
//
//     Eigen::Vector3d r_eci = sv.get_state().value().r;
//     Eigen::Vector3d v_eci = sv.get_state().value().v;
//     Eigen::Matrix3d dcm_eci2eqw = coordinateFrames::eci2eqw(r_eci, v_eci, mu);
//
//     // Compute total acceleration from force models
//     Eigen::Vector3d perturbedAcceleration_eci = computePerturbAcceleration(t, sv);
//     Eigen::Vector3d perturbedAcceleration_pqw = dcm_eci2eqw * perturbedAcceleration_eci;
//     double F_E = perturbedAcceleration_pqw(0);
//     double F_Q = perturbedAcceleration_pqw(1);
//     double F_W = perturbedAcceleration_pqw(2);
//
//
//     double X_constant = 1.0 + h_e * h_e + k_e * k_e;
//     double s_constant = std::sqrt(1.0 - a_f * a_f - a_g * a_g);
//     double c1 = 1.0 - ((a_f * a_f) / (1.0 + s_constant));
//     double c2 = (a_f * a_g) / (1.0 + s_constant);
//     double c3 = 1.0 - ((a_g * a_g) / (1.0 + s_constant));
//     double W_star = 1.0 - a_f * std::cos(eccLon) - a_g * std::sin(eccLon);
//     double s_constant_star = 1.0 / (1.0 + s_constant);
//
//     double dsma_dt = 2 * std::sqrt(sma * sma / (mu * W_star)) * (
//         ((c2 * std::cos(eccLon) - c1 * std::sin(eccLon)) * F_E
//         + (c1 * std::cos(eccLon) - c2 * std::sin(eccLon)) * F_Q)
//     );
//
//     double dhe_dt = std::sqrt(sma / mu) * (1 / s_constant) * (X_constant / 2) * (-a_f + c3 * std::cos(eccLon) + c2 * std::sin(eccLon)) * F_W;
//     double dke_dt = std::sqrt(sma / mu) * (1 / s_constant) * (X_constant / 2) * (-a_g + c2 * std::cos(eccLon) + c1 * std::sin(eccLon)) * F_W;
//
//     double daf_dt = std::sqrt(sma / mu) * (1.0 / W_star) * (
//         ((a_g * (c1 * std::cos(eccLon) - c2 * std::sin(eccLon))
//         + (-1.0 + a_f * a_f + 2.0 * c2 * c2) * std::cos(eccLon) * std::sin(eccLon)
//         + c3 * (s_constant + a_g * a_g * s_constant_star) * std::sin(eccLon) * std::sin(eccLon)
//         - c1 * c2 * std::cos(eccLon) * std::cos(eccLon)) * F_E
//         + (-a_f * (s_constant + c1) * std::cos(eccLon) + a_g * (-s_constant + a_f * a_f * s_constant_star) * std::sin(eccLon)
//         + c2 * (-1 - s_constant + 2 * c1) * std::cos(eccLon) * std::sin(eccLon)
//         + (s_constant - c2 * c2) * std::sin(eccLon) * std::sin(eccLon)
//         + (1.0 + 2.0 * s_constant - a_f * a_f - c1 * c1) * std::cos(eccLon) * std::cos(eccLon)) * F_Q)
//     ) + (2.0 * a_g / X_constant) * (k_e * dhe_dt - h_e * dke_dt);
//
//     double dag_dt = std::sqrt(sma / mu) * (1.0 / W_star) * (
//         ((a_f * (s_constant - a_g * a_g * s_constant_star) * std::cos(eccLon)
//         + a_g * (s_constant + c3) * std::sin(eccLon) + (-s_constant + c2 * c2) * std::cos(eccLon) * std::cos(eccLon)
//         + (-1.0 - 2.0 * s_constant + a_g * a_g + c3 * c3) * std::sin(eccLon) * std::sin(eccLon)
//         + c2 * (1.0 + s_constant - 2.0 * c3) * std::cos(eccLon) * std::sin(eccLon)) * F_E
//         + (a_f * (c2 * std::cos(eccLon) - c3 * std::sin(eccLon))
//         - c3 * (s_constant + a_f * a_f * s_constant_star) * std::cos(eccLon) * std::cos(eccLon)
//         + c2 * c3 * std::sin(eccLon) * std::sin(eccLon)
//         + (1.0 - a_g * a_g - 2.0 * c2 * c2) * std::cos(eccLon) * std::sin(eccLon)) * F_Q)
//     ) - (2.0 * a_f / X_constant) * (k_e * dhe_dt - h_e * dke_dt);
//
//     double deccLon_dt = std::sqrt(mu / (sma * sma * sma)) * (1.0 / W_star) * (
//     1.0 + (sma * sma / mu) * s_constant_star * (
//         (-s_constant * (1.0 + W_star) * std::cos(eccLon) + a_f * s_constant_star * (s_constant + (2.0 + s_constant) * W_star)
//         + a_g * s_constant_star * W_star * (a_f * std::sin(eccLon) - a_g * std::cos(eccLon))) * F_E
//         + (-s_constant * (1.0 + W_star) * std::sin(eccLon) + a_g * s_constant_star * (s_constant + (2 + s_constant) * W_star)
//         - a_f * s_constant_star * W_star * (a_f * std::sin(eccLon) - a_g * std::cos(eccLon))) * F_Q
//         + ((s_constant + 1.0) / s_constant) * W_star * (k_e * a_f - h_e * a_g
//         + (h_e * c2 - k_e * c3) * std::cos(eccLon)
//         + (h_e * c1 - k_e * c2) * std::sin(eccLon)) * F_W)
//     );
//
//     // Create and return the derivative of the state (dydt)
//     Vector6d dydt;
//     dydt << dsma_dt, daf_dt, dag_dt, dhe_dt, dke_dt, deccLon_dt;
//     return dydt;
// }
//
//
// void VariationOfParameters::updateState(SpaceVehicle& sv, double t) {
//     using Vector6d = Eigen::Matrix<double, 6, 1>;
//
//     SpaceVehicle::State svState = sv.get_state().value();
//     const double sma = svState.sma;
//     const double ecc = svState.ecc;
//     const double inc = svState.inc;
//     const double node = svState.node;
//     const double argP = svState.argP;
//     const double truA = svState.truA;
//
//     const double dt = computeStepSize(sma, ecc);
//     const double dt_max = computeMaxStep(sma, ecc);
//     const double dt_min = computeMinStep(sma, ecc);
//
//     const double tStart = sv.get_time().value();
//     const double tEnd = t;
//
//     Vector6d stateVector(6);
//
//     if (ecc < 0.001 || inc < 0.001) {
//         stateConversions::EquinoctialElements equin = stateConversions::Kep2Equinoctial(sma, ecc, inc, node, argP, truA);
//         stateVector << equin.sma, equin.a_f, equin.a_g, equin.h_e, equin.k_e, equin.eccLon;
//
//         integrator->integrate(
//         [this, &sv](const double t, const Vector6d& stateVector) -> Eigen::VectorXd {
//             return this->equationsOfMotion(t, stateVector, sv);
//             },
//         tStart, tEnd, stateVector, dt, dt_min, dt_max
//         );
//
//         stateConversions::KeplerianElements kep = stateConversions::Equinoctial2Kep(stateVector(0), stateVector(1), stateVector(2), stateVector(3), stateVector(4), stateVector(5));
//         stateConversions::CartesianState cart = stateConversions::Kep2Cart(kep.sma, kep.ecc, kep.inc, kep.node, kep.argP, kep.truA);
//         Eigen::Vector3d new_r = cart.r;
//         Eigen::Vector3d new_v = cart.v;
//
//         sv.set_rv(t, new_r, new_v);
//
//     } else {
//         stateVector << sma, ecc, inc, node, argP, truA;
//
//         integrator->integrate(
//         [this, &sv](const double t, const Vector6d& stateVector) -> Eigen::VectorXd {
//             return this->equationsOfMotion(t, stateVector, sv);
//             },
//     tStart, tEnd, stateVector, dt, dt_min, dt_max
//         );
//
//         stateConversions::CartesianState cart = stateConversions::Kep2Cart(stateVector(0), stateVector(1), stateVector(2), stateVector(3), stateVector(4), stateVector(5));
//         Eigen::Vector3d new_r = cart.r;
//         Eigen::Vector3d new_v = cart.v;
//
//         sv.set_rv(t, new_r, new_v);
//
//     }
//
// }


