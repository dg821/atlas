//
// Created by Douglas Garza on 1/14/25.
//

#include "../../../src/astro/basic_astrodynamics/stateConversions.h"
#include "../../../src/math/utilities/MathFunctions.h"
#include <iostream>
#include <cmath>
#include <gtest/gtest.h>
#include <Eigen/Dense>

class MeanTruEccAnomParameterizedTest : public ::testing::TestWithParam<std::tuple<double, double>> {

protected:
    const double angleTol = 1e-12;

    void SetUp() override {
        auto [eccAnom, ecc] = GetParam();

    }

};


INSTANTIATE_TEST_CASE_P(
    ElementVariations,
    MeanTruEccAnomParameterizedTest,
    ::testing::Combine(
        ::testing::Values(0.0, M_PI/4.0, M_PI/2.0, 3.0 * M_PI/4.0, M_PI, 5.0 * M_PI/4.0, 3.0*M_PI/2.0, 7.0 * M_PI/4.0),            // Eccentric Anomaly in rad
        ::testing::Values(0.0001, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 0.0)             // Eccentricity
    )
);

TEST_P(MeanTruEccAnomParameterizedTest, EccMeanTest) {
    double eccAnom = std::get<0>(GetParam());
    double ecc = std::get<1>(GetParam());

    double meanAnom = stateConversions::eccA2MeanA(eccAnom, ecc);
    double eccAnomBack = stateConversions::meanA2EccA(meanAnom, ecc);

    EXPECT_NEAR(eccAnomBack, eccAnom, angleTol);
}

TEST_P(MeanTruEccAnomParameterizedTest, EccTruTest) {
    double eccAnom = std::get<0>(GetParam());
    double ecc = std::get<1>(GetParam());

    double trueAnom = stateConversions::EccA2truA(eccAnom, ecc);
    double eccAnomBack = stateConversions::truA2EccA(trueAnom, ecc);

    EXPECT_NEAR(eccAnomBack, eccAnom, angleTol);
}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}