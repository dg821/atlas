//
// Created by Douglas Garza on 1/16/25.
//
#include "../../../src/astro/basic_astrodynamics/stateConversions.h"
#include "../../../src/math/utilities/mathFunctions.h"
#include "../../../src/math/UniversalConstants.h"
#include "../../../src/astro/coordinate_frames/coordinateFrames.h"
#include <iostream>
#include <cmath>
#include <gtest/gtest.h>
#include <Eigen/Dense>

class CoordinateSatFramesParameterizedTest : public ::testing::TestWithParam<std::tuple<double, double, double, double, double, double>> {

protected:
    const double tol = 1e-8;
    const double detTol = 1e-14;

    void SetUp() override {
        auto [sma, ecc, inc, node, argP, truA] = GetParam();
    }

};


INSTANTIATE_TEST_CASE_P(
    ElementVariations,
    CoordinateSatFramesParameterizedTest,
    ::testing::Combine(
        ::testing::Values(6700.0, 7000.0, 42164, 120000),
        ::testing::Values(0.001, 0.001, 0.1, 0.3, 0.5, 0.9),             // Eccentricity
        ::testing::Values(0.001, M_PI/4.0, M_PI/2.0 - 0.001),
        ::testing::Values(-M_PI/4.0, -M_PI/2.0, -3.0 * M_PI/4.0, -M_PI, -5.0 * M_PI/4.0, -3.0*M_PI/2.0, -7.0 * M_PI/4.0, 0.0, M_PI/4.0, M_PI/2.0, 3.0 * M_PI/4.0, M_PI, 5.0 * M_PI/4.0, 3.0*M_PI/2.0, 7.0 * M_PI/4.0),
        ::testing::Values(-M_PI/4.0, -M_PI/2.0, -3.0 * M_PI/4.0, -M_PI, -5.0 * M_PI/4.0, -3.0*M_PI/2.0, -7.0 * M_PI/4.0, 0.0, M_PI/4.0, M_PI/2.0, 3.0 * M_PI/4.0, M_PI, 5.0 * M_PI/4.0, 3.0*M_PI/2.0, 7.0 * M_PI/4.0),
        ::testing::Values(-M_PI/4.0, -M_PI/2.0, -3.0 * M_PI/4.0, -M_PI, -5.0 * M_PI/4.0, -3.0*M_PI/2.0, -7.0 * M_PI/4.0, 0.0, M_PI/4.0, M_PI/2.0, 3.0 * M_PI/4.0, M_PI, 5.0 * M_PI/4.0, 3.0*M_PI/2.0, 7.0 * M_PI/4.0)
    )
);

TEST_P(CoordinateSatFramesParameterizedTest, eci2rsw_test) {
    double sma = std::get<0>(GetParam());
    double ecc = std::get<1>(GetParam());
    double inc = std::get<2>(GetParam());
    double node = std::get<3>(GetParam());
    double argP = std::get<4>(GetParam());
    double truA = std::get<5>(GetParam());

    stateConversions::CartesianState cart = stateConversions::Kep2Cart( sma, ecc, inc, node, argP, truA);

    Eigen::Matrix3d dcm = coordinateFrames::eci2rsw(cart.r, cart.v);

    Eigen::Matrix3d inverse = dcm.inverse();
    Eigen::Matrix3d transpose = dcm.transpose();

    // Conditions for SO(3)
    EXPECT_NEAR(dcm.determinant(), 1, detTol);
    ASSERT_TRUE(inverse.isApprox(transpose));

    // Get back the same value
    Eigen::Vector3d doubleRot_r = dcm.transpose() * (dcm * cart.r);
    Eigen::Vector3d doubleRot_v = dcm.transpose() * (dcm * cart.v);

    ASSERT_TRUE(doubleRot_r.isApprox(cart.r));
    ASSERT_TRUE(doubleRot_v.isApprox(cart.v));

    Eigen::Vector3d h = cart.r.cross(cart.v);
    Eigen::Vector3d r_in_rsw = dcm * cart.r;
    Eigen::Vector3d h_in_rsw = dcm * h;
    Eigen::Vector3d w_in_rsw = dcm * (h.cross(cart.r));

    // vectors are correct
    ASSERT_TRUE((r_in_rsw.normalized()).isApprox(Eigen::Vector3d{1, 0, 0}));
    ASSERT_TRUE((h_in_rsw.normalized()).isApprox(Eigen::Vector3d{0, 0, 1}));
    ASSERT_TRUE((w_in_rsw.normalized()).isApprox(Eigen::Vector3d{0, 1, 0}));

}


TEST_P(CoordinateSatFramesParameterizedTest, eci2pqw_test) {
    double sma = std::get<0>(GetParam());
    double ecc = std::get<1>(GetParam());
    double inc = std::get<2>(GetParam());
    double node = std::get<3>(GetParam());
    double argP = std::get<4>(GetParam());
    double truA = std::get<5>(GetParam());

    double mu = UniversalConstants::EarthParams::MU;

    stateConversions::CartesianState cart = stateConversions::Kep2Cart( sma, ecc, inc, node, argP, truA);

    Eigen::Matrix3d dcm = coordinateFrames::eci2pqw(cart.r, cart.v);

    Eigen::Matrix3d inverse = dcm.inverse();
    Eigen::Matrix3d transpose = dcm.transpose();

    // Conditions for SO(3)
    EXPECT_NEAR(dcm.determinant(), 1, detTol);
    ASSERT_TRUE(inverse.isApprox(transpose));

    // Get back the same value
    Eigen::Vector3d doubleRot_r = dcm.transpose() * (dcm * cart.r);
    Eigen::Vector3d doubleRot_v = dcm.transpose() * (dcm * cart.v);

    ASSERT_TRUE(doubleRot_r.isApprox(cart.r));
    ASSERT_TRUE(doubleRot_v.isApprox(cart.v));

    Eigen::Vector3d eccVec = ((cart.v.norm() * cart.v.norm() - mu / cart.r.norm()) * cart.r - (cart.r.dot(cart.v)) * cart.v) / mu;

    Eigen::Vector3d eccVec_in_pqw = dcm * eccVec;
    Eigen::Vector3d h = cart.r.cross(cart.v);
    Eigen::Vector3d h_in_pqw = dcm * h;
    Eigen::Vector3d w_in_pqw = dcm * (h.cross(eccVec));

    Eigen::Vector3d eccVecUnit =  eccVec.normalized();

    // vectors are correct
    ASSERT_TRUE((eccVec_in_pqw.normalized()).isApprox(Eigen::Vector3d{1, 0, 0}));
    ASSERT_TRUE((h_in_pqw.normalized()).isApprox(Eigen::Vector3d{0, 0, 1}));
    ASSERT_TRUE((w_in_pqw.normalized()).isApprox(Eigen::Vector3d{0, 1, 0}));

}


TEST_P(CoordinateSatFramesParameterizedTest, eci2ntw_test) {
    double sma = std::get<0>(GetParam());
    double ecc = std::get<1>(GetParam());
    double inc = std::get<2>(GetParam());
    double node = std::get<3>(GetParam());
    double argP = std::get<4>(GetParam());
    double truA = std::get<5>(GetParam());

    double mu = UniversalConstants::EarthParams::MU;

    stateConversions::CartesianState cart = stateConversions::Kep2Cart( sma, ecc, inc, node, argP, truA);

    Eigen::Matrix3d dcm = coordinateFrames::eci2ntw(cart.r, cart.v);

    Eigen::Matrix3d inverse = dcm.inverse();
    Eigen::Matrix3d transpose = dcm.transpose();

    // Conditions for SO(3)
    EXPECT_NEAR(dcm.determinant(), 1, detTol);
    ASSERT_TRUE(inverse.isApprox(transpose));

    // Get back the same value
    Eigen::Vector3d doubleRot_r = dcm.transpose() * (dcm * cart.r);
    Eigen::Vector3d doubleRot_v = dcm.transpose() * (dcm * cart.v);

    ASSERT_TRUE(doubleRot_r.isApprox(cart.r));
    ASSERT_TRUE(doubleRot_v.isApprox(cart.v));

    Eigen::Vector3d v_in_rsw = dcm * cart.v;
    Eigen::Vector3d h = cart.r.cross(cart.v);
    Eigen::Vector3d h_in_rsw = dcm * h;
    Eigen::Vector3d n_in_rsw = dcm * (cart.v.cross(h));

    // vectors are correct
    ASSERT_TRUE((v_in_rsw.normalized()).isApprox(Eigen::Vector3d{1, 0, 0}));
    ASSERT_TRUE((h_in_rsw.normalized()).isApprox(Eigen::Vector3d{0, 1, 0}));
    ASSERT_TRUE((n_in_rsw.normalized()).isApprox(Eigen::Vector3d{0, 0, 1}));

}

TEST_P(CoordinateSatFramesParameterizedTest, eci2eqw_test) {
    double sma = std::get<0>(GetParam());
    double ecc = std::get<1>(GetParam());
    double inc = std::get<2>(GetParam());
    double node = std::get<3>(GetParam());
    double argP = std::get<4>(GetParam());
    double truA = std::get<5>(GetParam());

    stateConversions::CartesianState cart = stateConversions::Kep2Cart( sma, ecc, inc, node, argP, truA);

    Eigen::Matrix3d dcm = coordinateFrames::eci2eqw(cart.r, cart.v);

    Eigen::Matrix3d inverse = dcm.inverse();
    Eigen::Matrix3d transpose = dcm.transpose();

    // Conditions for SO(3)
    EXPECT_NEAR(dcm.determinant(), 1, detTol);
    ASSERT_TRUE(inverse.isApprox(transpose));

    // Get back the same value
    Eigen::Vector3d doubleRot_r = dcm.transpose() * (dcm * cart.r);
    Eigen::Vector3d doubleRot_v = dcm.transpose() * (dcm * cart.v);

    ASSERT_TRUE(doubleRot_r.isApprox(cart.r));
    ASSERT_TRUE(doubleRot_v.isApprox(cart.v));

}


TEST_P(CoordinateSatFramesParameterizedTest, pqw2rsw_test) {

    double truA = std::get<5>(GetParam());

    Eigen::Matrix3d dcm = coordinateFrames::pqw2rsw(truA);

    Eigen::Matrix3d inverse = dcm.inverse();
    Eigen::Matrix3d transpose = dcm.transpose();

    // Conditions for SO(3)
    EXPECT_NEAR(dcm.determinant(), 1, detTol);
    ASSERT_TRUE(inverse.isApprox(transpose));

}



int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}