set(PROPAGATORS_LIB_NAME propagators_lib)

set(propagators_HEADER_FILES
        OrbitPropagator.h
)

set(propagators_SOURCE_FILES
        OrbitPropagator.cpp
)

add_library(${PROPAGATORS_LIB_NAME}
        ${propagators_HEADER_FILES}
        ${propagators_SOURCE_FILES}
)

set_target_properties(${PROPAGATORS_LIB_NAME} PROPERTIES
        LINKER_LANGUAGE CXX
)

target_include_directories(${PROPAGATORS_LIB_NAME}
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(${PROPAGATORS_LIB_NAME}
        PUBLIC
        perturbation_methods_lib
        integrators_lib
        force_models_lib
)