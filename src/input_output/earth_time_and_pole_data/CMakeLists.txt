set(TIME_AND_POLE_DATA_LIB_NAME extract_time_and_pole_data_lib)

set(time_and_pole_data_HEADER_FILES
        ExtractTimeAndPoleData.h
        FundamentalArguments.h
        PrecessionNutationModel.h
        PrecessionNutation_YData.h
        PrecessionNutation_sData.h
        PrecessionNutation_XData.h
)

set(time_and_pole_data_SOURCE_FILES
        ExtractTimeAndPoleData.cpp
)

add_library(${TIME_AND_POLE_DATA_LIB_NAME}
        ${time_and_pole_data_HEADER_FILES}
        ${time_and_pole_data_SOURCE_FILES}
)

set_target_properties(${TIME_AND_POLE_DATA_LIB_NAME} PROPERTIES
        LINKER_LANGUAGE CXX
)

target_include_directories(${TIME_AND_POLE_DATA_LIB_NAME}
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)
