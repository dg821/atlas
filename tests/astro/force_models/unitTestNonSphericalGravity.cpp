//
// Created by Douglas Garza on 2/9/25.
//

//
// Created by Douglas Garza on 1/28/25.
//
#include "../../../src/astro/basic_astrodynamics/stateConversions.h"
#include "../../../src/math/utilities/mathFunctions.h"
#include "../../../src/math/UniversalConstants.h"
#include "../../../src/astro/propagators/OrbitPropagator.h"
#include "../../../src/astro/perturbation_methods/PerturbationMethods.h"
#include "../../../src/math/integrators/Dopri87.h"
#include "../../../src/astro/force_models/ForceModel.h"
#include <iostream>
#include <cmath>
#include <gtest/gtest.h>
#include <Eigen/Dense>


class NonSphericalGravityParameterizedTest : public ::testing::TestWithParam<std::tuple<double, double, double, double, double, double>> {

protected:
    const double tol = 1e-9;

    void SetUp() override {
        auto [sma, ecc, inc, node, argP, truA] = GetParam();

    }

};

INSTANTIATE_TEST_CASE_P(
    ElementVariations,
    NonSphericalGravityParameterizedTest,
    ::testing::Combine(
        ::testing::Values(6700.0, 7000.0, 7200.0),
        ::testing::Values(0.01, 0.1, 0.3),             // Eccentricity
        ::testing::Values(0.001, M_PI/4.0),
        ::testing::Values(0.0),
        ::testing::Values( 0.0),
        ::testing::Values( 0.0)
    )

    // Test cases contain a variety of orbit sizes, eccentricities (near circular and near parabolic), inclinations (near zero and near 90), node, argP, and truA
);

TEST_P(NonSphericalGravityParameterizedTest, J2TestSunSync) {
    double sma = std::get<0>(GetParam());
    double ecc = std::get<1>(GetParam());
    double inc = std::get<2>(GetParam());
    double node = std::get<3>(GetParam());
    double argP = std::get<4>(GetParam());
    double truA = std::get<5>(GetParam());

    double mu = UniversalConstants::EarthParams::MU;
    double j2 = UniversalConstants::EarthParams::J2;
    double rEq = UniversalConstants::EarthParams::RADIUS_EQ;

    double per = 2 * M_PI * std::sqrt(sma * sma * sma / mu);
    double t0 = 0;
    double t1 = t0 + per;

    double n = std::sqrt(mu / (sma * sma * sma));
    double nodeRateSunSync = -1.991063853e-7;            // rad/s
    double incSS = std::acos((-2.0 * std::pow(sma, 7.0/2.0) * nodeRateSunSync * ((1.0 - ecc*ecc) * (1.0 - ecc*ecc))) / (3.0 * rEq * rEq * j2 * std::sqrt(mu)));

    double p = sma * (1.0 - ecc * ecc);

    double apseRate = (3.0 * n * rEq * rEq * j2 / (4.0 * p * p)) * (4.0 - 5.0 * std::sin(incSS) * std::sin(incSS));

    SpaceVehicle sv = SpaceVehicle();
    sv.set_orb(t0, sma, ecc, incSS, node, argP, truA);

    OrbitPropagator propagator = OrbitPropagator();

    auto cowell = std::make_unique<Cowell>();
    propagator.setPerturbationMethod(std::move(cowell));

    auto dopri87 = std::make_unique<Dopri87>();
    propagator.setIntegrationMethod(std::move(dopri87));

    auto twoBody = std::make_unique<TwoBody>();
    auto nonSpherical_j2Only = std::make_unique<NonSphericalGravity>(2, 0);
    auto simpleJ2 = std::make_unique<J2SimpleTest>();
    propagator.addForceModel(std::move(twoBody));
    propagator.addForceModel(std::move(simpleJ2));

    propagator.propagate_to_time(sv, t1);

    double desNode2SS = 2*M_PI - nodeRateSunSync * (t1 - t0);
    double desApse2SS = 2*M_PI - apseRate * (t1 - t0);

    auto state = sv.get_state();

    double sma2 = state->sma;
    double ecc2 = state->ecc;
    double inc2 = state->inc;
    double node2 = state->node;
    double argP2 = state->argP;
    double truA2 = state->truA;

    EXPECT_NEAR(sma, sma2, tol);
    EXPECT_NEAR(ecc, ecc2, tol);
    EXPECT_NEAR(incSS, inc2, tol);
    EXPECT_NEAR(node2, desNode2SS, tol);
    EXPECT_NEAR(argP2, desApse2SS, tol);
    EXPECT_NEAR(truA, truA2, tol);

}

TEST_P(NonSphericalGravityParameterizedTest, J2SimpleVsJ2NonSpherical) {
    double sma = std::get<0>(GetParam());
    double ecc = std::get<1>(GetParam());
    double inc = std::get<2>(GetParam());
    double node = std::get<3>(GetParam());
    double argP = std::get<4>(GetParam());
    double truA = std::get<5>(GetParam());

    double mu = UniversalConstants::EarthParams::MU;
    double j2 = UniversalConstants::EarthParams::J2;
    double rEq = UniversalConstants::EarthParams::RADIUS_EQ;

    double per = 2 * M_PI * std::sqrt(sma * sma * sma / mu);
    double t0 = 0;
    double t1 = t0 + per;

    double n = std::sqrt(mu / (sma * sma * sma));
    double nodeRateSunSync = -1.991063853e-7;            // rad/s
    double incSS = std::acos((-2.0 * std::pow(sma, 7.0/2.0) * nodeRateSunSync * ((1.0 - ecc*ecc) * (1.0 - ecc*ecc))) / (3.0 * rEq * rEq * j2 * std::sqrt(mu)));

    double p = sma * (1.0 - ecc * ecc);

    double apseRate = (3.0 * n * rEq * rEq * j2 / (4.0 * p * p)) * (4.0 - 5.0 * std::sin(incSS) * std::sin(incSS));

    SpaceVehicle sv = SpaceVehicle();
    sv.set_orb(t0, sma, ecc, incSS, node, argP, truA);

    OrbitPropagator propagator = OrbitPropagator();

    auto cowell = std::make_unique<Cowell>();
    propagator.setPerturbationMethod(std::move(cowell));

    auto dopri87 = std::make_unique<Dopri87>();
    propagator.setIntegrationMethod(std::move(dopri87));

    auto twoBody = std::make_unique<TwoBody>();
    auto nonSpherical_j2Only = std::make_unique<NonSphericalGravity>(2, 0);
    auto simpleJ2 = std::make_unique<J2SimpleTest>();
    propagator.addForceModel(std::move(twoBody));
    propagator.addForceModel(std::move(simpleJ2));

    propagator.propagate_to_time(sv, t1);

    auto state = sv.get_state();

    double sma2Simple = state->sma;
    double ecc2Simple = state->ecc;
    double inc2Simple = state->inc;
    double node2Simple = state->node;
    double argP2Simple = state->argP;
    double truA2Simple = state->truA;

    sv.set_orb(t0, sma, ecc, incSS, node, argP, truA);
    propagator.removeForceModel("SimpleJ2");
    propagator.addForceModel(std::move(nonSpherical_j2Only));

    propagator.propagate_to_time(sv, t1);

    auto state2 = sv.get_state();

    double smaJ2 = state2->sma;
    double eccJ2 = state2->ecc;
    double incJ2  = state2->inc;
    double nodeJ2  = state2->node;
    double argPJ2  = state2->argP;
    double truAJ2  = state2->truA;


    EXPECT_NEAR(sma2Simple, smaJ2, tol);
    EXPECT_NEAR(ecc2Simple, eccJ2, tol);
    EXPECT_NEAR(inc2Simple, incJ2, tol);
    EXPECT_NEAR(node2Simple, nodeJ2, tol);
    EXPECT_NEAR(argP2Simple, argPJ2, tol);
    EXPECT_NEAR(truA2Simple, truAJ2, tol);

}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}