set(PLANETARY_EPHEMERIDES_LIB_NAME planetary_ephemerides_lib)

set(planetary_ephemerides_HEADER_FILES
        meeusEphemeris.h
)

add_library(${PLANETARY_EPHEMERIDES_LIB_NAME} INTERFACE)

target_include_directories(${PLANETARY_EPHEMERIDES_LIB_NAME}
        INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}
)