//
// Created by Douglas Garza on 1/22/25.
//
#include "../../../src/astro/geodetic_model/geodeticModel.h"
#include "../../../src/astro/basic_astrodynamics/stateConversions.h"
#include "../../../src/math/utilities/mathFunctions.h"
#include "../../../src/math/UniversalConstants.h"
#include "../../../src/astro/coordinate_frames/coordinateFrames.h"
#include <iostream>
#include <cmath>
#include <gtest/gtest.h>
#include <Eigen/Dense>


class GeodeticModelParameterizedTest : public ::testing::TestWithParam<std::tuple<Eigen::Vector3d, Eigen::Vector3d>> {

protected:
    const double tol = 1e-14;

    void SetUp() override {
        auto [r, v] = GetParam();
    }

};

INSTANTIATE_TEST_CASE_P(
    DateVariations,
    GeodeticModelParameterizedTest,
        ::testing::Values(std::make_tuple(Eigen::Vector3d {7000, 0, 0}, Eigen::Vector3d{0, 9, 0 }))
);

TEST_P(GeodeticModelParameterizedTest, single_case_itrf_to_tirs) {

    // // Vallado test
    Eigen::Vector3d r_tirs_truth = {-1033.4750312, 7901.3055856, 6380.3445327};

    Eigen::Vector3d r = {-1033.4793830, 7901.2952754, 6380.3565958};
    Eigen::Vector3d v = {-3.225636520, -2.872451450, 5.531924446};

    int year = 2004; int month = 4; int day = 6; int hour = 7; int minute = 51; double second = 28.386009;
    timeConversions::GregorianDate date {year, month, day, hour, minute, second};

    double timej2k = timeConversions::date_to_timej2k(date.year, date.month, date.day, date.hour, date.minute, date.second);

    Eigen::Matrix3d dcm_itrf2tirs = geodeticModel::getDCM_itrf2tirs(timej2k);

    Eigen::Vector3d r_tirs = dcm_itrf2tirs * r;

    ASSERT_TRUE(r_tirs.isApprox(r_tirs_truth));

}

TEST_P(GeodeticModelParameterizedTest, single_case_tirs_to_cirs) {

    // // Vallado test
    Eigen::Vector3d r_cirs_truth = {5100.0184047, 6122.7863648, 6380.3445327};

    Eigen::Vector3d r = {-1033.4750312, 7901.3055856, 6380.3445327};
    // Eigen::Vector3d v = {-3.225636520, -2.872451450, 5.531924446};

    int year = 2004; int month = 4; int day = 6; int hour = 7; int minute = 51; double second = 28.386009;
    timeConversions::GregorianDate date {year, month, day, hour, minute, second};

    double timej2k = timeConversions::date_to_timej2k(date.year, date.month, date.day, date.hour, date.minute, date.second);

    Eigen::Matrix3d dcm_tirs2cirs = geodeticModel::getDCM_tirs2cirs(timej2k);

    Eigen::Vector3d r_cirs = dcm_tirs2cirs * r;

    ASSERT_TRUE(r_cirs.isApprox(r_cirs_truth));

}

TEST_P(GeodeticModelParameterizedTest, single_case_cirs_to_gcrf) {

    Eigen::Vector3d r_gcrf_truth = {5102.508959, 6123.011403, 6378.136925};

    Eigen::Vector3d r = {5100.0184047, 6122.7863648, 6380.3445327};

    int year = 2004; int month = 4; int day = 6; int hour = 7; int minute = 51; double second = 28.386009;
    timeConversions::GregorianDate date {year, month, day, hour, minute, second};

    double timej2k = timeConversions::date_to_timej2k(date.year, date.month, date.day, date.hour, date.minute, date.second);

    Eigen::Matrix3d dcm_cirs2gcrf = geodeticModel::getDCM_cirs2gcrf_lofi(timej2k);

    Eigen::Vector3d r_gcrf = dcm_cirs2gcrf * r;

    ASSERT_TRUE(r_gcrf.isApprox(r_gcrf_truth));

}

TEST_P(GeodeticModelParameterizedTest, single_case_ecef_to_eci) {

    Eigen::Vector3d r_gcrf_truth = {5102.508959, 6123.011403, 6378.136925};
    Eigen::Vector3d v_gcrf_truth = {-4.74322016, 0.79053650, 5.53375573};

    Eigen::Vector3d r_ecef = {-1033.4793830, 7901.2952754, 6380.3565958};
    Eigen::Vector3d v_ecef = {-3.225636520, -2.872451450, 5.531924446};

    int year = 2004; int month = 4; int day = 6; int hour = 7; int minute = 51; double second = 28.386009;
    timeConversions::GregorianDate date {year, month, day, hour, minute, second};

    double timej2k = timeConversions::date_to_timej2k(date.year, date.month, date.day, date.hour, date.minute, date.second);

    auto res = geodeticModel::ecef2eci(r_ecef, v_ecef, timej2k);
    Eigen::Vector3d r_eci = std::get<0>(res);
    Eigen::Vector3d v_eci = std::get<1>(res);

    ASSERT_TRUE(r_eci.isApprox(r_gcrf_truth));
    ASSERT_TRUE(v_eci.isApprox(v_gcrf_truth));

}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}