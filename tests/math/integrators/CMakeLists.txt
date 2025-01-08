add_executable(unitTestDopri87 unitTestDopri87.cpp)

# Link libraries to the test target
target_link_libraries(unitTestDopri87
        PRIVATE
        GTest::gtest_main
        Eigen3::Eigen
        integrators_lib
        math_lib
        utilities_lib
        basic_astrodynamics_lib
)




# Set include directories for the test target
target_include_directories(unitTestDopri87
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../../..
)




add_test(NAME unitTestDopri87 COMMAND unitTestDopri87)
