add_executable(unitTestAngleConversions unitTestAngleConversions.cpp)
add_executable(unitTestOrbitStateConversions unitTestOrbitStateConversions.cpp)

# Link libraries to the test target
target_link_libraries(unitTestAngleConversions
        PRIVATE
        GTest::gtest_main
        Eigen3::Eigen
        utilities_lib
        integrators_lib
        math_lib
        basic_astrodynamics_lib
)

# Link libraries to the test target
target_link_libraries(unitTestOrbitStateConversions
        PRIVATE
        GTest::gtest_main
        Eigen3::Eigen
        utilities_lib
        integrators_lib
        math_lib
        basic_astrodynamics_lib
)




# Set include directories for the test target
target_include_directories(unitTestAngleConversions
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../../..
)

# Set include directories for the test target
target_include_directories(unitTestOrbitStateConversions
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../../..
)