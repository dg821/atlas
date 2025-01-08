set(SPACEVEHICLE_LIB_NAME space_vehicle_lib)

set(SpaceVehicle_HEADER_FILES
        SpaceVehicle.h
)

set(SpaceVehicle_SOURCE_FILES
        SpaceVehicle.cpp
)

add_library(space_vehicle_lib
        ${SpaceVehicle_HEADER_FILES}
        ${SpaceVehicle_SOURCE_FILES}
)

set_target_properties(${SPACEVEHICLE_LIB_NAME} PROPERTIES
        LINKER_LANGUAGE CXX
)

target_include_directories(${SPACEVEHICLE_LIB_NAME}
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)
