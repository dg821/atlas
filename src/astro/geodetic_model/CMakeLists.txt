set(GEODETIC_MODEL_LIB_NAME geodetic_model_lib)

set(geodetic_model_HEADER_FILES
        GeodeticModel.h
)

set(geodetic_model_SOURCE_FILES
        GeodeticModel.cpp
)

add_library(${GEODETIC_MODEL_LIB_NAME}
        ${geodetic_model_HEADER_FILES}
        ${geodetic_model_SOURCE_FILES}
)

set_target_properties(${GEODETIC_MODEL_LIB_NAME} PROPERTIES
        LINKER_LANGUAGE CXX
)

target_include_directories(${GEODETIC_MODEL_LIB_NAME}
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)