set(PERTURBATION_METHODS_LIB_NAME perturbation_methods_lib)

set(perturbation_methods_HEADER_FILES
        PerturbationMethods.h
)

set(perturbation_methods_SOURCE_FILES
        PerturbationMethods.cpp
)

add_library(${PERTURBATION_METHODS_LIB_NAME}
        ${perturbation_methods_HEADER_FILES}
        ${perturbation_methods_SOURCE_FILES}
)

set_target_properties(${PERTURBATION_METHODS_LIB_NAME} PROPERTIES
        LINKER_LANGUAGE CXX
)

target_include_directories(${PERTURBATION_METHODS_LIB_NAME}
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(${PERTURBATION_METHODS_LIB_NAME}
        PUBLIC
        space_vehicle_lib
        basic_astrodynamics_lib
)



