//
// Created by Douglas Garza on 2/2/25.
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


class NonsphericalGravityClassParameterizedTest : public ::testing::Test {

protected:
    const double tol = 1e-9;


};


TEST_F(NonsphericalGravityClassParameterizedTest, ClassTest) {

    NonSphericalGravity();

}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}