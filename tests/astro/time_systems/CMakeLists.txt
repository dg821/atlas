set(TEST_TARGETS
        unitTestTimeConversions
)

foreach(TARGET ${TEST_TARGETS})
    add_executable(${TARGET} ${TARGET}.cpp)

    target_link_libraries(${TARGET}
            PRIVATE
            GTest::gtest_main
            Eigen3::Eigen
            coordinate_frames_lib
            utilities_lib
            math_lib
            basic_astrodynamics_lib
            time_conversions_lib
            extract_time_and_pole_data_lib
    )

    target_include_directories(${TARGET}
            PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/../../..
    )
endforeach()