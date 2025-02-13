set(FORCE_MODELS_LIB_NAME force_models_lib)

set(force_models_HEADER_FILES
        ForceModel.h
        ExponentialDragLookupTable.h
)

set(force_models_SOURCE_FILES
        ForceModel.cpp
)

add_library(${FORCE_MODELS_LIB_NAME}
        ${force_models_HEADER_FILES}
        ${force_models_SOURCE_FILES}
)

set_target_properties(${FORCE_MODELS_LIB_NAME} PROPERTIES
        LINKER_LANGUAGE CXX
)

target_include_directories(${FORCE_MODELS_LIB_NAME}
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(${FORCE_MODELS_LIB_NAME}
        PUBLIC
        space_vehicle_lib
)


