//
// Created by Douglas Garza on 1/14/25.
//

#include "../../../src/astro/basic_astrodynamics/stateConversions.h"
#include "../../../src/math/utilities/MathFunctions.h"
#include "../../../src/math/UniversalConstants.h"
#include <iostream>
#include <cmath>
#include <gtest/gtest.h>
#include <Eigen/Dense>


class StateConversionsParameterizedTest : public ::testing::TestWithParam<std::tuple<double, double, double, double, double, double>> {

protected:
    const double angleTol = 1e-7;
    const double eccTol = 1e-8;
    const double smaTol = 1e-8;

    void SetUp() override {
        auto [sma, ecc, inc, node, argP, truA] = GetParam();

    }

};


INSTANTIATE_TEST_CASE_P(
    ElementVariations,
    StateConversionsParameterizedTest,
    ::testing::Combine(
        ::testing::Values(6700.0, 7000.0, 42164, 120000),
        ::testing::Values(0.001, 0.001, 0.1, 0.3, 0.5, 0.9),             // Eccentricity
        ::testing::Values(0.001, M_PI/4.0, M_PI/2.0 - 0.001),
        ::testing::Values(0.0, M_PI/4.0, M_PI/2.0, 3.0*M_PI/4.0, M_PI, 5.0*M_PI/4.0, 3.0*M_PI/2.0, 7.0*M_PI/4.0),
        ::testing::Values(0.0, M_PI/4.0, M_PI/2.0, 3.0*M_PI/4.0, M_PI, 5.0*M_PI/4.0, 3.0*M_PI/2.0, 7.0*M_PI/4.0),
        ::testing::Values(0.0, M_PI/4.0, M_PI/2.0, 3.0*M_PI/4.0, M_PI, 5.0*M_PI/4.0, 3.0*M_PI/2.0, 7.0*M_PI/4.0)
    )

    // Test cases contain a variety of orbit sizes, eccentricities (near circular and near parabolic), inclinations (near zero and near 90), node, argP, and truA
);

TEST_P(StateConversionsParameterizedTest, CartesianKeplerianTest) {

    // Test: convert keplerian to cartesian then back to keplerian and compare output and input

    double sma = std::get<0>(GetParam());
    double ecc = std::get<1>(GetParam());
    double inc = std::get<2>(GetParam());
    double node = std::get<3>(GetParam());
    double argP = std::get<4>(GetParam());
    double truA = std::get<5>(GetParam());

    double mu = UniversalConstants::EarthParams::MU;

    stateConversions::CartesianState cart = stateConversions::Kep2Cart(sma, ecc, inc, node, argP, truA, mu);
    stateConversions::KeplerianElements kep = stateConversions::Cart2Kep(cart.r, cart.v, mu);

    // Ignore scenarios that compare 0 and 2pi
    if (abs(kep.node - 2.0*M_PI) < angleTol) {
        kep.node = 0.0;
    }
    if (abs(kep.argP - 2.0*M_PI) < angleTol) {
        kep.argP = 0.0;
    }
    if (abs(kep.truA - 2.0*M_PI) < angleTol) {
        kep.truA = 0.0;
    }

    EXPECT_NEAR(sma, kep.sma, smaTol);
    EXPECT_NEAR(ecc, kep.ecc, eccTol);
    EXPECT_NEAR(inc, kep.inc, angleTol);
    EXPECT_NEAR(node, kep.node, angleTol);
    EXPECT_NEAR(argP, kep.argP, angleTol);
    EXPECT_NEAR(truA, kep.truA, angleTol);
}

TEST_P(StateConversionsParameterizedTest, EqunioctalConversionTest) {
    // Test: convert keplerian to equinoctial then back to keplerian and compare output and input


    double sma = std::get<0>(GetParam());
    double ecc = std::get<1>(GetParam());
    double inc = std::get<2>(GetParam());
    double node = std::get<3>(GetParam());
    double argP = std::get<4>(GetParam());
    double truA = std::get<5>(GetParam());


    stateConversions::EquinoctialElements equin = stateConversions::Kep2Equinoctial(sma,  ecc,  inc,  node,  argP,  truA);
    double sma2 = equin.sma;
    double a_f = equin.a_f;
    double a_g = equin.a_g;
    double h_e = equin.h_e;
    double k_e = equin.k_e;
    double eccLon = equin.eccLon;

    stateConversions::KeplerianElements kep = stateConversions::Equinoctial2Kep( sma2,  a_f,  a_g,  h_e,  k_e,  eccLon);

    // Ignore scenarios that compare 0 and 2pi

    if (abs(kep.node - 2.0*M_PI) < angleTol) {
        kep.node = 0.0;
    }
    if (abs(kep.argP - 2.0*M_PI) < angleTol) {
        kep.argP = 0.0;
    }
    if (abs(kep.truA - 2.0*M_PI) < angleTol) {
        kep.truA = 0.0;
    }

    EXPECT_NEAR(sma, kep.sma, smaTol);
    EXPECT_NEAR(ecc, kep.ecc, eccTol);
    EXPECT_NEAR(inc, kep.inc, angleTol);
    EXPECT_NEAR(node, kep.node, angleTol);
    EXPECT_NEAR(argP, kep.argP, angleTol);
    EXPECT_NEAR(truA, kep.truA, angleTol);

}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}