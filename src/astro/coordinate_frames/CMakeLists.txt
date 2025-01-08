set(COORDINATE_FRAMES_LIB_NAME coordinate_frames_lib)

set(coordinate_frames_HEADER_FILES
        CoordinateFrames.h
)

set(coordinate_frames_SOURCE_FILES
        CoordinateFrames.cpp
)

add_library(${COORDINATE_FRAMES_LIB_NAME}
        ${coordinate_frames_HEADER_FILES}
        ${coordinate_frames_SOURCE_FILES}
)

set_target_properties(${COORDINATE_FRAMES_LIB_NAME} PROPERTIES
        LINKER_LANGUAGE CXX
)

target_include_directories(${COORDINATE_FRAMES_LIB_NAME}
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)