#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include "../../math/utilities/NumericalDiff.h"

class DifferentiatorTest : public ::testing::Test {
protected:
    // Test points from -2π to 2π
    std::vector<double> test_points;
    static constexpr double tolerance = 1e-6;

    void SetUp() override {
        const int num_points = 20;
        for (int i = 0; i < num_points; ++i) {
            test_points.push_back(-2.0 * M_PI + 4.0 * M_PI * i / (num_points - 1));
        }
    }

    // Helper to test a function and its known derivative
    void TestFunction(
        const std::function<double(double)>& f,
        const std::function<double(double)>& f_prime,
        const std::string& func_name
    ) {
        for (double x : test_points) {
            double error_estimate;
            double numerical_deriv = NumericalDiff::differentiate(f, x, tolerance);
            double exact_deriv = f_prime(x);
            double abs_error = std::abs(numerical_deriv - exact_deriv);

            EXPECT_LE(abs_error, tolerance)
                << "Failed for " << func_name << " at x = " << x << "\n"
                << "Numerical derivative: " << numerical_deriv << "\n"
                << "Exact derivative: " << exact_deriv << "\n"
                << "Error estimate: " << error_estimate;
        }
    }
};

// Test sin(x) -> cos(x)
TEST_F(DifferentiatorTest, SineFunction) {
    auto f = [](double x) { return std::sin(x); };
    auto f_prime = [](double x) { return std::cos(x); };
    TestFunction(f, f_prime, "sin(x)");
}

// Test cos(x) -> -sin(x)
TEST_F(DifferentiatorTest, CosineFunction) {
    auto f = [](double x) { return std::cos(x); };
    auto f_prime = [](double x) { return -std::sin(x); };
    TestFunction(f, f_prime, "cos(x)");
}

// Test exp(x) -> exp(x)
TEST_F(DifferentiatorTest, ExponentialFunction) {
    auto f = [](double x) { return std::exp(x); };
    auto f_prime = [](double x) { return std::exp(x); };
    TestFunction(f, f_prime, "exp(x)");
}

// Test x^3 -> 3x^2
TEST_F(DifferentiatorTest, CubicFunction) {
    auto f = [](double x) { return x * x * x; };
    auto f_prime = [](double x) { return 3 * x * x; };
    TestFunction(f, f_prime, "x^3");
}

// Test 1/x -> -1/x^2
TEST_F(DifferentiatorTest, ReciprocalFunction) {
    auto f = [](double x) { return 1.0 / x; };
    auto f_prime = [](double x) { return -1.0 / (x * x); };

    // Skip x = 0 for this test
    std::vector<double> nonzero_points;
    std::copy_if(test_points.begin(), test_points.end(),
                 std::back_inserter(nonzero_points),
                 [](double x) { return std::abs(x) > 0.1; });

    for (double x : nonzero_points) {
        double error_estimate;
        double numerical_deriv = NumericalDiff::differentiate(f, x, tolerance);
        double exact_deriv = f_prime(x);
        double abs_error = std::abs(numerical_deriv - exact_deriv);

        EXPECT_LE(abs_error, tolerance)
            << "Failed for 1/x at x = " << x << "\n"
            << "Numerical derivative: " << numerical_deriv << "\n"
            << "Exact derivative: " << exact_deriv << "\n"
            << "Error estimate: " << error_estimate;
    }
}

// Test error estimation accuracy
// TEST_F(DifferentiatorTest, ErrorEstimation) {
//     auto f = [](double x) { return std::sin(x); };
//     auto f_prime = [](double x) { return std::cos(x); };
//
//     for (double x : test_points) {
//         double error_estimate;
//         double numerical_deriv = NumericalDiff::differentiate(f, x, error_estimate, tolerance);
//         double exact_deriv = f_prime(x);
//         double actual_error = std::abs(numerical_deriv - exact_deriv);
//
//         // Error estimate should be within an order of magnitude
//         EXPECT_LE(actual_error, 10 * error_estimate)
//             << "Error estimate too small at x = " << x;
//         EXPECT_LE(error_estimate, 100 * actual_error)
//             << "Error estimate too large at x = " << x;
//     }
// }

// Test pathological function
TEST_F(DifferentiatorTest, HighlyOscillatoryFunction) {
    // f(x) = sin(1/x) - highly oscillatory near x = 0
    auto f = [](double x) { return std::sin(1.0 / x); };

    std::vector<double> safe_points = {1.0, 2.0, 3.0, 4.0, 5.0};
    for (double x : safe_points) {
        double error_estimate;
        EXPECT_NO_THROW({
            NumericalDiff::differentiate(f, x, tolerance);
        }) << "Failed to handle oscillatory function at x = " << x;
    }
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}