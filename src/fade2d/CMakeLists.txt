cmake_minimum_required(VERSION 3.10)
project(fade2d)

# Determine the Visual Studio version being used
set(VS_VERSION "v141")

# Define the imported library target for the debug configuration
add_library(fade2d_debug STATIC IMPORTED)
set_target_properties(fade2d_debug PROPERTIES
    IMPORTED_LOCATION "${CMAKE_CURRENT_SOURCE_DIR}/x64/fade2D_x64_${VS_VERSION}_Debug.lib"
    INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/include_fade2d"	
)

# Define the imported library target for the release configuration
add_library(fade2d_release STATIC IMPORTED)
set_target_properties(fade2d_release PROPERTIES
    IMPORTED_LOCATION "${CMAKE_CURRENT_SOURCE_DIR}/x64/fade2D_x64_${VS_VERSION}_Release.lib"
    INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/include_fade2d"
)

# Create an alias target for convenience
add_library(fade2d INTERFACE)
target_link_libraries(fade2d INTERFACE
    $<$<CONFIG:Debug>:fade2d_debug>
    $<$<CONFIG:Release>:fade2d_release>
)
