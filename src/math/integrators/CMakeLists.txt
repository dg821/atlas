set(INTEGRATORS_LIB_NAME integrators_lib)

set(integrators_HEADER_FILES
        NumericalIntegrator.h
        Dopri87.h
)

set(integrators_SOURCE_FILES
        NumericalIntegrator.cpp
        Dopri87.cpp
)

add_library(${INTEGRATORS_LIB_NAME}
        ${integrators_HEADER_FILES}
        ${integrators_SOURCE_FILES}
)

set_target_properties(${INTEGRATORS_LIB_NAME} PROPERTIES
        LINKER_LANGUAGE CXX
)

target_include_directories(${INTEGRATORS_LIB_NAME}
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

