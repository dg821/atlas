set(TIME_CONVERSIONS_LIB_NAME time_conversions_lib)

set(time_conversions_HEADER_FILES
        TimeConversions.h
)

set(time_conversions_SOURCE_FILES
        TimeConversions.cpp
)

add_library(${TIME_CONVERSIONS_LIB_NAME}
        ${time_conversions_HEADER_FILES}
        ${time_conversions_SOURCE_FILES}
)

set_target_properties(${TIME_CONVERSIONS_LIB_NAME} PROPERTIES
        LINKER_LANGUAGE CXX
)

target_include_directories(${TIME_CONVERSIONS_LIB_NAME}
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)
