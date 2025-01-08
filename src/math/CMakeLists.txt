add_subdirectory("integrators")
add_subdirectory("utilities")

set(MATH_LIB_NAME math_lib)

set(math_HEADER_FILES
        UniversalConstants.h
)

# Create an interface target for the header-only library
add_library(${MATH_LIB_NAME} INTERFACE)

# Add include directories for this library
target_include_directories(${MATH_LIB_NAME}
        INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# If you want to specify usage of any particular language version (optional)
set_target_properties(${MATH_LIB_NAME} PROPERTIES
        LINKER_LANGUAGE CXX
)