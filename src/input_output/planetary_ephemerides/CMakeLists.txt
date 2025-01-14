set(PLANETARY_EPHEMERIDES_LIB_NAME planetary_ephemerides_lib)

set(planetary_ephemerides_HEADER_FILES
        meeusEphemeris.h
)

add_library(${PLANETARY_EPHEMERIDES_LIB_NAME}
        ${planetary_ephemerides_HEADER_FILES}
)

set_target_properties(${PLANETARY_EPHEMERIDES_LIB_NAME} PROPERTIES
        LINKER_LANGUAGE CXX
)

target_include_directories(${PLANETARY_EPHEMERIDES_LIB_NAME}
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)
