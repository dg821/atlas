set(TEST_TARGETS
        unitTestNumericalDiff
)

foreach(TARGET ${TEST_TARGETS})
    add_executable(${TARGET} ${TARGET}.cpp)

    target_link_libraries(${TARGET}
            PRIVATE
            GTest::gtest_main
            Eigen3::Eigen
            utilities_lib
            math_lib
    )

    target_include_directories(${TARGET}
            PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/../../..
    )
endforeach()