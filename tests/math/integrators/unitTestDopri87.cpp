#include "../../../src/math/integrators/Dopri87.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <functional>
#include <fstream>
#include <stdexcept>
#include <gtest/gtest.h>
#include "../../../src/astro/propagators/OrbitPropagator.h"
#include "../../../src/astro/basic_astrodynamics/stateConversions.h"
#include "../../../src/math/UniversalConstants.h"

// Forward declarations
double calculateEnergy(const NumericalIntegrator::StateVector& state);
void WriteTestResultsToCSV(const std::string& filename,
                           const std::vector<std::pair<double, NumericalIntegrator::StateVector>>& solution,
                           const std::tuple<double, double, double, double, double, double>& params,
                           const std::chrono::duration<double>& execution_time);
void CalculateAndRecordErrorMetrics(const std::string& filename,
                                    const NumericalIntegrator::StateVector& initial_state,
                                    const NumericalIntegrator::StateVector& final_state,
                                    const std::tuple<double, double, double, double, double, double>& params);
double calculateInitialStepSize(double sma, double ecc, double mu);
double getMinStepSize(double sma, double ecc, double mu);
double getMaxStepSize(double sma, double ecc, double mu);

// 3D orbital motion equations
class Orbit3D {
public:
    static NumericalIntegrator::StateVector f(double t, const NumericalIntegrator::StateVector& y) {
        const double mu = UniversalConstants::EarthParams::MU;

        double rMag = std::sqrt(y(0)*y(0) + y(1)*y(1) + y(2)*y(2));
        double r3 = std::pow(rMag, 3);

        NumericalIntegrator::StateVector dydx(6);
        dydx(0) = y(3);
        dydx(1) = y(4);
        dydx(2) = y(5);
        dydx(3) = -mu * y(0) / r3;
        dydx(4) = -mu * y(1) / r3;
        dydx(5) = -mu * y(2) / r3;

        return dydx;
    }
};

// Helper function to calculate total energy of the system
double calculateEnergy(const NumericalIntegrator::StateVector& state) {
    const double mu = UniversalConstants::EarthParams::MU;

    Eigen::Vector3d vNow = {state(3), state(4), state(5)};
    Eigen::Vector3d rNow = {state(0), state(1), state(2)};

    double specE = vNow.squaredNorm() / 2 - mu / rNow.norm();

    return specE;
}

void WriteTestResultsToCSV(const std::string& filename,
                           const std::vector<std::pair<double, NumericalIntegrator::StateVector>>& solution,
                           const std::tuple<double, double, double, double, double, double>& params,
                           const std::chrono::duration<double>& execution_time) {
    std::ofstream outFile(filename, std::ios::app);
    if (!outFile.is_open()) {
        throw std::runtime_error("Unable to open file: " + filename);
    }

    auto [sma, ecc, inc, node, argP, truA] = params;

    outFile.seekp(0, std::ios::end);
    if (outFile.tellp() == 0) {
        outFile << "sma,ecc,inc,node,argP,truA,execution_time,time,x,y,z,vx,vy,vz\n";
    }

    for (const auto& [t, state] : solution) {
        outFile << sma << "," << ecc << "," << inc << "," << node << "," << argP << "," << truA << ","
                << execution_time.count() << "," << t << ","
                << state(0) << "," << state(1) << "," << state(2) << ","
                << state(3) << "," << state(4) << "," << state(5) << "\n";
    }

    outFile.close();
}

void CalculateAndRecordErrorMetrics(const std::string& filename,
                                    const NumericalIntegrator::StateVector& initial_state,
                                    const NumericalIntegrator::StateVector& final_state,
                                    const std::tuple<double, double, double, double, double, double>& params) {
    std::ofstream outFile(filename, std::ios::app);
    if (!outFile.is_open()) {
        throw std::runtime_error("Unable to open file: " + filename);
    }

    auto [sma, ecc, inc, node, argP, truA] = params;

    double pos_error = (final_state.head(3) - initial_state.head(3)).norm();
    double vel_error = (final_state.tail(3) - initial_state.tail(3)).norm();

    outFile.seekp(0, std::ios::end);
    if (outFile.tellp() == 0) {
        outFile << "sma,ecc,inc,node,argP,truA,position_error,velocity_error\n";
    }

    outFile << sma << "," << ecc << "," << inc << "," << node << "," << argP << "," << truA << ","
            << pos_error << "," << vel_error << "\n";

    outFile.close();
}

double calculateInitialStepSize(double sma, double ecc, double mu) {
    double period = 2.0 * M_PI * std::sqrt(std::pow(sma, 3) / mu);
    double fraction = std::pow(10, log10(sma) - 1);
    return (period / fraction) * std::pow(1.0 - ecc, 2);
}

double getMinStepSize(double sma, double ecc, double mu) {
    double period = 2.0 * M_PI * std::sqrt(std::pow(sma, 3) / mu);
    double fraction = std::pow(10, log10(sma) + 2);
    return (period / fraction) * std::pow(1.0 - ecc, 3);
}

double getMaxStepSize(double sma, double ecc, double mu) {
    double period = 2.0 * M_PI * std::sqrt(std::pow(sma, 3) / mu);
    double fraction = std::pow(10, floor(log10(sma) - 1));
    return (period / fraction) * std::pow(1.0 - ecc, 2);
}

class Dopri87ParameterizedTest : public ::testing::TestWithParam<std::tuple<double, double, double, double, double, double>> {
protected:
    Dopri87 integrator;
    NumericalIntegrator::StateVector initial_state;
    double per;
    const double mu = UniversalConstants::EarthParams::MU;
    const double Req = UniversalConstants::EarthParams::RADIUS_EQ;
    double h, h_min, h_max;

    void SetUp() override {
        auto [sma, ecc, inc, node, argP, truA] = GetParam();

        stateConversions::CartesianState cart = stateConversions::Kep2Cart(sma, ecc, inc, node, argP, truA, mu);
        initial_state = NumericalIntegrator::StateVector(6);
        initial_state << cart.r(0), cart.r(1), cart.r(2), cart.v(0), cart.v(1), cart.v(2);

        std::cout << "SetUp: initial_state size = " << initial_state.size() << std::endl;

        per = 2 * M_PI * std::sqrt(std::pow(sma, 3) / mu);

        h = calculateInitialStepSize(sma, ecc, mu);
        h_min = getMinStepSize(sma, ecc, mu);
        h_max = getMaxStepSize(sma, ecc, mu);
    }
};

INSTANTIATE_TEST_CASE_P(
    OrbitVariations,
    Dopri87ParameterizedTest,
    ::testing::Combine(
        ::testing::Values(7000.0, 26050.0, 42164.0, 120000.0),  // sma in km
        ::testing::Values(0.001, 0.3),         // ecc
        ::testing::Values(0.01, M_PI/4, M_PI/2 - 0.1),          // inc
        ::testing::Values(0.0, -M_PI/4),     // node
        ::testing::Values(0.0, -M_PI/4),     // argP
        ::testing::Values(0.0, M_PI)      // truA
    )
);

TEST(Dopri87ParameterizedTest, ConstructorTest) {
    Dopri87 integrator;
    EXPECT_EQ(integrator.getName(), "Dormand-Prince 8(7)");
    EXPECT_EQ(integrator.getOrder(), 8);
}

TEST_P(Dopri87ParameterizedTest, SingleStepTest) {
    double t = 0.0;
    integrator.setTolerances(1e-12, 1e-12);

    // Perform a single step
    NumericalIntegrator::StateVector state_new = integrator.step(Orbit3D::f, t, initial_state, h);

    // Check that the new state vector has the same size as the initial state
    ASSERT_EQ(state_new.size(), initial_state.size()) << "State vector sizes do not match";

    // Check that the new state is different from the initial state
    EXPECT_FALSE(state_new.isApprox(initial_state, 1e-10)) << "New state is too close to initial state";

}

TEST_P(Dopri87ParameterizedTest, IntegrationTest) {
    integrator.setTolerances(1e-12, 1e-12);
    double t_start = 0.0;
    double t_end = per/2.0;

    auto start_time = std::chrono::high_resolution_clock::now();
    auto solution = integrator.integrate(Orbit3D::f, t_start, t_end, initial_state, h, h_min, h_max);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto execution_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);

    WriteTestResultsToCSV("integration_test_results_2.csv", solution, GetParam(), execution_time);

    auto result = solution.back();
    CalculateAndRecordErrorMetrics("integration_test_errors_2.csv", initial_state, result.second, GetParam());

    const double tolerance = 1e-7;
    for (int i = 0; i < 6; ++i) {
        EXPECT_NEAR(result.second(i), initial_state(i), tolerance);
    }
}

TEST_P(Dopri87ParameterizedTest, EnergyConservationTest) {
    double t_start = 0.0;
    double t_end = 0.5 * per;

    auto start_time = std::chrono::high_resolution_clock::now();
    auto solution = integrator.integrate(Orbit3D::f, t_start, t_end, initial_state, h, h_min, h_max);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto execution_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);

    double initial_energy = calculateEnergy(initial_state);

    const double tolerance = 1e-10;
    for (const auto& [t, state] : solution) {
        double energy = calculateEnergy(state);
        EXPECT_NEAR(energy, initial_energy, tolerance);
    }

    WriteTestResultsToCSV("energy_conservation_test_2.csv", solution, GetParam(), execution_time);

    std::cout << "Energy Conservation Test execution time: " << execution_time.count() << " seconds" << std::endl;
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}