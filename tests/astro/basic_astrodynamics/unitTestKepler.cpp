//
// Created by Douglas Garza on 1/16/25.
//

#include "../../../src/astro/basic_astrodynamics/stateConversions.h"
#include "../../../src/math/utilities/MathFunctions.h"
#include "../../../src/math/UniversalConstants.h"
#include "../../../src/astro/basic_astrodynamics/Kepler.h"
#include <iostream>
#include <cmath>
#include <gtest/gtest.h>
#include <Eigen/Dense>

class KeplerProblemParameterizedTest : public ::testing::TestWithParam<std::tuple<double, double, double, double, double, double>> {

protected:
    const double tol = 1e-6;
    const double tolHighEcc = 1e-5;

    void SetUp() override {
        auto [sma, ecc, inc, node, argP, truA] = GetParam();
    }

};

INSTANTIATE_TEST_CASE_P(
    ElementVariations,
    KeplerProblemParameterizedTest,
    ::testing::Combine(
    ::testing::Values(6700.0, 7000.0, 42164, 120000),
    ::testing::Values(0.001, 0.001, 0.1, 0.3, 0.5, 0.9),             // Eccentricity
    ::testing::Values(0.001, M_PI/4.0, M_PI/2.0 - 0.001),
    ::testing::Values(-M_PI/4.0, -M_PI/2.0, -3.0 * M_PI/4.0, -M_PI, -5.0 * M_PI/4.0, -3.0*M_PI/2.0, -7.0 * M_PI/4.0, 0.0, M_PI/4.0, M_PI/2.0, 3.0 * M_PI/4.0, M_PI, 5.0 * M_PI/4.0, 3.0*M_PI/2.0, 7.0 * M_PI/4.0),
    ::testing::Values(-M_PI/4.0, -M_PI/2.0, -3.0 * M_PI/4.0, -M_PI, -5.0 * M_PI/4.0, -3.0*M_PI/2.0, -7.0 * M_PI/4.0, 0.0, M_PI/4.0, M_PI/2.0, 3.0 * M_PI/4.0, M_PI, 5.0 * M_PI/4.0, 3.0*M_PI/2.0, 7.0 * M_PI/4.0),
    ::testing::Values(-M_PI/4.0, -M_PI/2.0, -3.0 * M_PI/4.0, -M_PI, -5.0 * M_PI/4.0, -3.0*M_PI/2.0, -7.0 * M_PI/4.0, 0.0, M_PI/4.0, M_PI/2.0, 3.0 * M_PI/4.0, M_PI, 5.0 * M_PI/4.0, 3.0*M_PI/2.0, 7.0 * M_PI/4.0)
    )
);



TEST_P(KeplerProblemParameterizedTest, ReturnAfterFullPeriodTest) {
    double sma = std::get<0>(GetParam());
    double ecc = std::get<1>(GetParam());
    double inc = std::get<2>(GetParam());
    double node = std::get<3>(GetParam());
    double argP = std::get<4>(GetParam());
    double truA = std::get<5>(GetParam());

    double mu = UniversalConstants::EarthParams::MU;

    stateConversions::CartesianState cart = stateConversions::Kep2Cart(sma, ecc, inc, node, argP, truA);

    Eigen::Vector3d r = cart.r;
    Eigen::Vector3d v = cart.v;


    double per = 2 * M_PI * std::sqrt(sma * sma * sma / mu);

    auto res = Kepler::solveKepler(r, v, per);

    Eigen::Vector3d rOut = res.first;
    Eigen::Vector3d vOut = res.second;

    if (ecc > 0.7) {
        EXPECT_NEAR(rOut(0), r(0), tol);
        EXPECT_NEAR(rOut(1), r(1), tol);
        EXPECT_NEAR(rOut(2), r(2), tol);
        EXPECT_NEAR(vOut(0), v(0), tol);
        EXPECT_NEAR(vOut(1), v(1), tol);
        EXPECT_NEAR(vOut(2), v(2), tol);
    } else {
        EXPECT_NEAR(rOut(0), r(0), tolHighEcc);
        EXPECT_NEAR(rOut(1), r(1), tolHighEcc);
        EXPECT_NEAR(rOut(2), r(2), tolHighEcc);
        EXPECT_NEAR(vOut(0), v(0), tolHighEcc);
        EXPECT_NEAR(vOut(1), v(1), tolHighEcc);
        EXPECT_NEAR(vOut(2), v(2), tolHighEcc);
    }



}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}