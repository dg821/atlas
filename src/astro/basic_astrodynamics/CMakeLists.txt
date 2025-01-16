set(BASIC_ASTRODYNAMICS_LIB_NAME basic_astrodynamics_lib)

set(basic_astrodynamics_HEADER_FILES
        stateConversions.h
        Kepler.h
)

set(basic_astrodynamics_SOURCE_FILES
        stateConversions.cpp
        Kepler.cpp
)

add_library(${BASIC_ASTRODYNAMICS_LIB_NAME}
        ${basic_astrodynamics_HEADER_FILES}
        ${basic_astrodynamics_SOURCE_FILES}
)

set_target_properties(${BASIC_ASTRODYNAMICS_LIB_NAME} PROPERTIES
        LINKER_LANGUAGE CXX
)

target_include_directories(${BASIC_ASTRODYNAMICS_LIB_NAME}
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(${BASIC_ASTRODYNAMICS_LIB_NAME}
        PUBLIC
        math_lib
        utilities_lib
        Eigen3::Eigen
)

message(STATUS "basic_astrodynamics_lib will be built in: ${CMAKE_CURRENT_BINARY_DIR}")
message(STATUS "basic_astrodynamics_lib source files: ${basic_astrodynamics_SOURCE_FILES}")