//
// Created by Douglas Garza on 2/12/25.
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


void writeOrbitDataToCSV(const std::vector<double>& timeVec, const std::vector<Eigen::Vector3d>& positionVec, const std::string& filename) {
    // Create full path for the output file
    std::filesystem::path currentPath = std::filesystem::current_path();

    std::filesystem::path outputDir = currentPath / "test_output";
    std::filesystem::create_directories(outputDir); // Create directory if it doesn't exist
    std::filesystem::path fullPath = outputDir / filename;

    std::cout << "Writing file to: " << fullPath << std::endl;

    // Open file with error checking
    std::ofstream outFile(fullPath);
    if (!outFile.is_open()) {
        throw std::runtime_error("Could not open file: " + fullPath.string());
    }

    // Write header
    outFile << "Time(s),X(km),Y(km),Z(km)\n";

    // Write data
    for (size_t i = 0; i < timeVec.size(); i++) {
        outFile << std::fixed << std::setprecision(6)  // Set precision for floating point
                << timeVec[i] << ","
                << positionVec[i](0) << ","
                << positionVec[i](1) << ","
                << positionVec[i](2) << "\n";
    }

    outFile.close();

    // Verify file was written
    if (!std::filesystem::exists(fullPath)) {
        throw std::runtime_error("Failed to write file: " + fullPath.string());
    }
}

// Test fixture class
class ForceModelStrengthTest : public ::testing::Test {
protected:
    // SetUp is called before each test
    void SetUp() override {
        // Initialize any test resources here
    }

    // TearDown is called after each test
    void TearDown() override {
        // Clean up any test resources here
    }
};

// TEST_F(ForceModelStrengthTest, Model_TwoBody) {
//     // setup tVec
//     double t0 = 0;
//     double tf = t0 + 5760 * UniversalConstants::SECONDS_PER_MINUTE;
//     double dt = 100;
//     int idx = 0;
//     double t = t0;
//     int size = static_cast<int>((tf - t0) / dt) + 1;
//     std::vector<double> tVec(size);
//     while (t <= tf) {
//         tVec[idx] = t;
//         t += dt;
//         idx++;
//     }
//
//     // setup orbit
//     double rp = UniversalConstants::EarthParams::RADIUS_EQ + 380;
//     double ra = UniversalConstants::EarthParams::RADIUS_EQ + 390;
//     double inc = 51.6 * UniversalConstants::D2R;
//     double sma = (rp + ra) / 2;
//     double ecc = (ra - rp) / (ra + rp);
//     double node = 0.0;
//     double argP = 0.0;
//     double truA = 0.0;
//
//     SpaceVehicle svTest = SpaceVehicle();
//     svTest.set_orb(t0, sma, ecc, inc, node, argP, truA);
//
//     OrbitPropagator testProp = OrbitPropagator();
//
//     auto cowell = std::make_unique<Cowell>();
//     testProp.setPerturbationMethod(std::move(cowell));
//
//     auto dopri87 = std::make_unique<Dopri87>();
//     testProp.setIntegrationMethod(std::move(dopri87));
//
//     auto twoBody = std::make_unique<TwoBody>();
//     auto nonSpherical_j2Only = std::make_unique<NonSphericalGravity>(2, 0);
//     testProp.addForceModel(std::move(twoBody));
//
//     std::vector<Eigen::Vector3d> testR_vec(tVec.size());
//     for (int i = 0; i < tVec.size(); i++) {
//         testProp.propagate_to_time(svTest, tVec[i]);
//
//         auto testState = svTest.get_state();
//
//         auto testR = testState->r;
//         testR_vec[i] = testR;
//
//         // std::cout << "t: " << tVec[i] << "       rVec: " << testR(0) << "  " << testR(1) << "  " << testR(2) << "  " << std::endl;
//     }
//
//     writeOrbitDataToCSV(tVec, testR_vec, "/Users/douglasgarza/CLionProjects/atlas/test_data/force_model_magnitude/two_body_pos_data.csv");
//
// }
//
// TEST_F(ForceModelStrengthTest, Model_2x0) {
//     // setup tVec
//     double t0 = 0;
//     double tf = t0 + 5760 * UniversalConstants::SECONDS_PER_MINUTE;
//     double dt = 100;
//     int idx = 0;
//     double t = t0;
//     int size = static_cast<int>((tf - t0) / dt) + 1;
//     std::vector<double> tVec(size);
//     while (t <= tf) {
//         tVec[idx] = t;
//         t += dt;
//         idx++;
//     }
//
//     // setup orbit
//     double rp = UniversalConstants::EarthParams::RADIUS_EQ + 380;
//     double ra = UniversalConstants::EarthParams::RADIUS_EQ + 390;
//     double inc = 51.6 * UniversalConstants::D2R;
//     double sma = (rp + ra) / 2;
//     double ecc = (ra - rp) / (ra + rp);
//     double node = 0.0;
//     double argP = 0.0;
//     double truA = 0.0;
//
//     SpaceVehicle svTest = SpaceVehicle();
//     svTest.set_orb(t0, sma, ecc, inc, node, argP, truA);
//
//     OrbitPropagator testProp = OrbitPropagator();
//
//     auto cowell = std::make_unique<Cowell>();
//     testProp.setPerturbationMethod(std::move(cowell));
//
//     auto dopri87 = std::make_unique<Dopri87>();
//     testProp.setIntegrationMethod(std::move(dopri87));
//
//     auto twoBody = std::make_unique<TwoBody>();
//     auto nonSpherical_j2Only = std::make_unique<NonSphericalGravity>(2, 0);
//     testProp.addForceModel(std::move(twoBody));
//     testProp.addForceModel(std::move(nonSpherical_j2Only));
//
//     std::vector<Eigen::Vector3d> testR_vec(tVec.size());;
//     for (int i = 0; i < tVec.size(); i++) {
//         testProp.propagate_to_time(svTest, tVec[i]);
//
//         auto testState = svTest.get_state();
//
//         auto testR = testState->r;
//         testR_vec[i] = testR;
//
//
//         // std::cout << "t: " << tVec[i] << "       rVec: " << testR(0) << "  " << testR(1) << "  " << testR(2) << "  " << std::endl;
//     }
//
//     writeOrbitDataToCSV(tVec, testR_vec, "/Users/douglasgarza/CLionProjects/atlas/test_data/force_model_magnitude/2x0_pos_data.csv");
//
//
// }
//

TEST_F(ForceModelStrengthTest, Model_4x0) {
    // setup tVec
    double t0 = 0;
    double tf = t0 + 5760 * UniversalConstants::SECONDS_PER_MINUTE;
    double dt = 100;
    int idx = 0;
    double t = t0;
    int size = static_cast<int>((tf - t0) / dt) + 1;
    std::vector<double> tVec(size);
    while (t <= tf) {
        tVec[idx] = t;
        t += dt;
        idx++;
    }

    // setup orbit
    double rp = UniversalConstants::EarthParams::RADIUS_EQ + 380;
    double ra = UniversalConstants::EarthParams::RADIUS_EQ + 390;
    double inc = 51.6 * UniversalConstants::D2R;
    double sma = (rp + ra) / 2;
    double ecc = (ra - rp) / (ra + rp);
    double node = 0.0;
    double argP = 0.0;
    double truA = 0.0;

    SpaceVehicle svTest = SpaceVehicle();
    svTest.set_orb(t0, sma, ecc, inc, node, argP, truA);

    OrbitPropagator testProp = OrbitPropagator();

    auto cowell = std::make_unique<Cowell>();
    testProp.setPerturbationMethod(std::move(cowell));

    auto dopri87 = std::make_unique<Dopri87>();
    testProp.setIntegrationMethod(std::move(dopri87));

    auto twoBody = std::make_unique<TwoBody>();
    auto nonSpherical_j23 = std::make_unique<NonSphericalGravity>(4, 0);
    testProp.addForceModel(std::move(twoBody));
    testProp.addForceModel(std::move(nonSpherical_j23));

    std::vector<Eigen::Vector3d> testR_vec(tVec.size());;
    for (int i = 0; i < tVec.size(); i++) {
        testProp.propagate_to_time(svTest, tVec[i]);

        auto testState = svTest.get_state();

        auto testR = testState->r;
        auto testV = testState->v;
        testR_vec[i] = testR;


        std::cout << "t: " << tVec[i] << "       rVec: " << testR(0) << "  " << testR(1) << "  " << testR(2) << ",        "
        << "vVec: " << testV(0) << "  " << testV(1) << "  " << testV(2) << "  " << std::endl;
    }

    writeOrbitDataToCSV(tVec, testR_vec, "/Users/douglasgarza/CLionProjects/atlas/test_data/force_model_magnitude/4x0_pos_data.csv");


}

// TEST_F(ForceModelStrengthTest, Model_Drag) {
//     // setup tVec
//     double t0 = 0;
//     double tf = t0 + 5760 * UniversalConstants::SECONDS_PER_MINUTE;
//     double dt = 100;
//     int idx = 0;
//     double t = t0;
//     int size = static_cast<int>((tf - t0) / dt) + 1;
//     std::vector<double> tVec(size);
//     while (t <= tf) {
//         tVec[idx] = t;
//         t += dt;
//         idx++;
//     }
//
//     // setup orbit
//     double rp = UniversalConstants::EarthParams::RADIUS_EQ + 380;
//     double ra = UniversalConstants::EarthParams::RADIUS_EQ + 390;
//     double inc = 51.6 * UniversalConstants::D2R;
//     double sma = (rp + ra) / 2;
//     double ecc = (ra - rp) / (ra + rp);
//     double node = 0.0;
//     double argP = 0.0;
//     double truA = 0.0;
//
//     SpaceVehicle svTest = SpaceVehicle();
//     svTest.set_orb(t0, sma, ecc, inc, node, argP, truA);
//
//     OrbitPropagator testProp = OrbitPropagator();
//
//     auto cowell = std::make_unique<Cowell>();
//     testProp.setPerturbationMethod(std::move(cowell));
//
//     auto dopri87 = std::make_unique<Dopri87>();
//     testProp.setIntegrationMethod(std::move(dopri87));
//
//     auto twoBody = std::make_unique<TwoBody>();
//     auto expoDrag = std::make_unique<DragExpForce>();
//     testProp.addForceModel(std::move(twoBody));
//     testProp.addForceModel(std::move(expoDrag));
//
//     std::vector<Eigen::Vector3d> testR_vec(tVec.size());;
//     for (int i = 0; i < tVec.size(); i++) {
//         testProp.propagate_to_time(svTest, tVec[i]);
//
//         auto testState = svTest.get_state();
//
//         auto testR = testState->r;
//         testR_vec[i] = testR;
//
//         std::cout << "t: " << tVec[i] << "       rVec: " << testR(0) << "  " << testR(1) << "  " << testR(2) << "  " << std::endl;
//     }
//
//     writeOrbitDataToCSV(tVec, testR_vec, "/Users/douglasgarza/CLionProjects/atlas/test_data/force_model_magnitude/drag_pos_data.csv");
//
// }

// TEST_F(ForceModelStrengthTest, Model_SRP) {
//     // setup tVec
//     double t0 = 0;
//     double tf = t0 + 5760 * UniversalConstants::SECONDS_PER_MINUTE;
//     double dt = 100;
//     int idx = 0;
//     double t = t0;
//     int size = static_cast<int>((tf - t0) / dt) + 1;
//     std::vector<double> tVec(size);
//     while (t <= tf) {
//         tVec[idx] = t;
//         t += dt;
//         idx++;
//     }
//
//     // setup orbit
//     double rp = UniversalConstants::EarthParams::RADIUS_EQ + 380;
//     double ra = UniversalConstants::EarthParams::RADIUS_EQ + 390;
//     double inc = 51.6 * UniversalConstants::D2R;
//     double sma = (rp + ra) / 2;
//     double ecc = (ra - rp) / (ra + rp);
//     double node = 0.0;
//     double argP = 0.0;
//     double truA = 0.0;
//
//     SpaceVehicle svTest = SpaceVehicle();
//     svTest.set_orb(t0, sma, ecc, inc, node, argP, truA);
//
//     OrbitPropagator testProp = OrbitPropagator();
//
//     auto cowell = std::make_unique<Cowell>();
//     testProp.setPerturbationMethod(std::move(cowell));
//
//     auto dopri87 = std::make_unique<Dopri87>();
//     testProp.setIntegrationMethod(std::move(dopri87));
//
//     auto twoBody = std::make_unique<TwoBody>();
//     auto srp = std::make_unique<SimpleSRPForce>();
//     testProp.addForceModel(std::move(twoBody));
//     testProp.addForceModel(std::move(srp));
//
//     std::vector<Eigen::Vector3d> testR_vec(tVec.size());;
//     for (int i = 0; i < tVec.size(); i++) {
//         testProp.propagate_to_time(svTest, tVec[i]);
//
//         auto testState = svTest.get_state();
//
//         auto testR = testState->r;
//         testR_vec[i] = testR;
//
//         std::cout << "t: " << tVec[i] << "       rVec: " << testR(0) << "  " << testR(1) << "  " << testR(2) << "  " << std::endl;
//     }
//
//     writeOrbitDataToCSV(tVec, testR_vec, "/Users/douglasgarza/CLionProjects/atlas/test_data/force_model_magnitude/srp_pos_data.csv");
//
// }

// TEST_F(ForceModelStrengthTest, Model_ThirdBody) {
//     // setup tVec
//     double t0 = 0;
//     double tf = t0 + 5760 * UniversalConstants::SECONDS_PER_MINUTE;
//     double dt = 100;
//     int idx = 0;
//     double t = t0;
//     int size = static_cast<int>((tf - t0) / dt) + 1;
//     std::vector<double> tVec(size);
//     while (t <= tf) {
//         tVec[idx] = t;
//         t += dt;
//         idx++;
//     }
//
//     // setup orbit
//     double rp = UniversalConstants::EarthParams::RADIUS_EQ + 380;
//     double ra = UniversalConstants::EarthParams::RADIUS_EQ + 390;
//     double inc = 51.6 * UniversalConstants::D2R;
//     double sma = (rp + ra) / 2;
//     double ecc = (ra - rp) / (ra + rp);
//     double node = 0.0;
//     double argP = 0.0;
//     double truA = 0.0;
//
//     SpaceVehicle svTest = SpaceVehicle();
//     svTest.set_orb(t0, sma, ecc, inc, node, argP, truA);
//
//     OrbitPropagator testProp = OrbitPropagator();
//
//     auto cowell = std::make_unique<Cowell>();
//     testProp.setPerturbationMethod(std::move(cowell));
//
//     auto dopri87 = std::make_unique<Dopri87>();
//     testProp.setIntegrationMethod(std::move(dopri87));
//
//     auto twoBody = std::make_unique<TwoBody>();
//     auto thirdBody = std::make_unique<LunisolarThirdBodyForce>();
//     testProp.addForceModel(std::move(twoBody));
//     testProp.addForceModel(std::move(thirdBody));
//
//     std::vector<Eigen::Vector3d> testR_vec(tVec.size());;
//     for (int i = 0; i < tVec.size(); i++) {
//         testProp.propagate_to_time(svTest, tVec[i]);
//
//         auto testState = svTest.get_state();
//
//         auto testR = testState->r;
//         testR_vec[i] = testR;
//
//         std::cout << "t: " << tVec[i] << "       rVec: " << testR(0) << "  " << testR(1) << "  " << testR(2) << "  " << std::endl;
//     }
//
//     writeOrbitDataToCSV(tVec, testR_vec, "/Users/douglasgarza/CLionProjects/atlas/test_data/force_model_magnitude/thirdBody_pos_data.csv");
//
// }

// TEST_F(ForceModelStrengthTest, Model_12x12) {
//     // setup tVec
//     double t0 = 0;
//     double tf = t0 + 5760 * UniversalConstants::SECONDS_PER_MINUTE;
//     double dt = 100;
//     int idx = 0;
//     double t = t0;
//     int size = static_cast<int>((tf - t0) / dt) + 1;
//     std::vector<double> tVec(size);
//     while (t <= tf) {
//         tVec[idx] = t;
//         t += dt;
//         idx++;
//     }
//
//     // setup orbit
//     double rp = UniversalConstants::EarthParams::RADIUS_EQ + 380;
//     double ra = UniversalConstants::EarthParams::RADIUS_EQ + 390;
//     double inc = 51.6 * UniversalConstants::D2R;
//     double sma = (rp + ra) / 2;
//     double ecc = (ra - rp) / (ra + rp);
//     double node = 0.0;
//     double argP = 0.0;
//     double truA = 0.0;
//
//     SpaceVehicle svTest = SpaceVehicle();
//     svTest.set_orb(t0, sma, ecc, inc, node, argP, truA);
//
//     OrbitPropagator testProp = OrbitPropagator();
//
//     auto cowell = std::make_unique<Cowell>();
//     testProp.setPerturbationMethod(std::move(cowell));
//
//     auto dopri87 = std::make_unique<Dopri87>();
//     testProp.setIntegrationMethod(std::move(dopri87));
//
//     auto twoBody = std::make_unique<TwoBody>();
//     auto nonSpherical_12x12 = std::make_unique<NonSphericalGravity>(5, 0);
//     testProp.addForceModel(std::move(twoBody));
//     testProp.addForceModel(std::move(nonSpherical_12x12));
//
//     std::vector<double> diffRMag;
//     for (int i = 0; i < tVec.size(); i++) {
//         testProp.propagate_to_time(svTest, tVec[i]);
//
//         auto testState = svTest.get_state();
//
//         auto testR = testState->r;
//
//         std::cout << "t: " << tVec[i] << "       rVec: " << testR(0) << "  " << testR(1) << "  " << testR(2) << "  " << std::endl;
//     }
//
// }

// TEST_F(ForceModelStrengthTest, Model_20x20) {
//     // setup tVec
//     double t0 = 0;
//     double tf = t0 + 5760 * UniversalConstants::SECONDS_PER_MINUTE;
//     double dt = 100;
//     int idx = 0;
//     double t = t0;
//     int size = static_cast<int>((tf - t0) / dt) + 1;
//     std::vector<double> tVec(size);
//     while (t <= tf) {
//         tVec[idx] = t;
//         t += dt;
//         idx++;
//     }
//
//     // setup orbit
//     double rp = UniversalConstants::EarthParams::RADIUS_EQ + 380;
//     double ra = UniversalConstants::EarthParams::RADIUS_EQ + 390;
//     double inc = 51.6 * UniversalConstants::D2R;
//     double sma = (rp + ra) / 2;
//     double ecc = (ra - rp) / (ra + rp);
//     double node = 0.0;
//     double argP = 0.0;
//     double truA = 0.0;
//
//     SpaceVehicle svTest = SpaceVehicle();
//     svTest.set_orb(t0, sma, ecc, inc, node, argP, truA);
//
//     OrbitPropagator testProp = OrbitPropagator();
//
//     auto cowell = std::make_unique<Cowell>();
//     testProp.setPerturbationMethod(std::move(cowell));
//
//     auto dopri87 = std::make_unique<Dopri87>();
//     testProp.setIntegrationMethod(std::move(dopri87));
//
//     auto twoBody = std::make_unique<TwoBody>();
//     auto nonSpherical_20x20 = std::make_unique<NonSphericalGravity>(20, 20);
//     testProp.addForceModel(std::move(twoBody));
//     testProp.addForceModel(std::move(nonSpherical_20x20));
//
//     std::vector<double> diffRMag;
//     for (int i = 0; i < tVec.size(); i++) {
//         testProp.propagate_to_time(svTest, tVec[i]);
//
//         auto testState = svTest.get_state();
//
//         auto testR = testState->r;
//
//         std::cout << "t: " << tVec[i] << "       rVec: " << testR(0) << "  " << testR(1) << "  " << testR(2) << "  " << std::endl;
//     }
//
// }


int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}