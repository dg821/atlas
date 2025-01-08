cmake_minimum_required(VERSION 3.29)
project(atlas)

set(CMAKE_CXX_STANDARD 20)


include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(/usr/local/include/eigen/3.4.0_1/include/eigen3)
include_directories(/Users/douglasgarza/googletest/googletest/include/gtest)

set(googletest_DIR "/Users/douglasgarza/googletest/googletest/cmake")
set(Eigen3_DIR "/usr/local/include/eigen/3.4.0_1/share/eigen3/cmake")

find_package(GTest REQUIRED)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_subdirectory(src/math)
add_subdirectory(src/astro)
add_subdirectory(tests/math)
add_subdirectory(src/input_ouput)


add_executable(atlas main.cpp
        src/astro/geodetic_model/GeodeticModel.h)

# Link libraries
target_link_libraries(atlas
        PRIVATE
        GTest::gtest
        Eigen3::Eigen
        math_lib
        basic_astrodynamics_lib
        space_vehicle_lib
        utilities_lib
        propagators_lib
        force_models_lib
        perturbation_methods_lib
        integrators_lib
        time_conversions_lib
        geodetic_model_lib
        coordinate_frames_lib
        extract_time_and_pole_data_lib
)

# Add include directories for the main target
target_include_directories(atlas
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

enable_testing()

