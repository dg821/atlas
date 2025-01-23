set(UTILITIES_LIB_NAME utilities_lib)

set(utilities_HEADER_FILES
    RotationMatrices.h
        MathFunctions.h
)

set(utilities_SOURCE_FILES
    RotationMatrices.cpp
        MathFunctions.cpp
)

add_library(${UTILITIES_LIB_NAME}
        ${utilities_HEADER_FILES}
        ${utilities_SOURCE_FILES}
)

set_target_properties(${UTILITIES_LIB_NAME} PROPERTIES
        LINKER_LANGUAGE CXX
)

target_include_directories(${UTILITIES_LIB_NAME}
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)