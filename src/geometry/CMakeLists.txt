cmake_minimum_required(VERSION 3.14)
Project(Geometry)

set(CMAKE_CONFIGURATION_TYPES Debug Release)

#add_compile_options(/bigobj) # for some Eigen solvers

# collect only the top-level Geometry files
file(GLOB Geometry_Src
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
)
file(GLOB Geometry_Includes
    "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
)

# collect only the Poisson folder files
file(GLOB Poisson_Files
    "${CMAKE_CURRENT_SOURCE_DIR}/Poisson/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/Poisson/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/Poisson/*.inl"
)

# build the library from both sets
add_library(${PROJECT_NAME}
    ${Geometry_Src}
    ${Geometry_Includes}
    ${Poisson_Files}
)

# Create VS filters
if (MSVC)
  source_group("Poisson" FILES ${Poisson_Files})
endif()

# Includes (so #include "Poisson/Octree.h" works)
target_include_directories(${PROJECT_NAME}
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    "${CMAKE_CURRENT_SOURCE_DIR}/../fade2d/include_fade2d"
)

# Link Geometry library with other libraries, including Fade2D
target_link_libraries(${PROJECT_NAME} PUBLIC pmp Utils fade2d poly2tri sawhney_mat triangle tet)

# Copy Fade2D DLLs to the output directory
set(VS_VERSION "v141")
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${CMAKE_CURRENT_SOURCE_DIR}/../fade2d/x64/fade2D_x64_${VS_VERSION}_Debug.dll"
    "$<TARGET_FILE_DIR:${PROJECT_NAME}>/fade2D_x64_${VS_VERSION}_Debug.dll"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${CMAKE_CURRENT_SOURCE_DIR}/../fade2d/x64/fade2D_x64_${VS_VERSION}_Release.dll"
    "$<TARGET_FILE_DIR:${PROJECT_NAME}>/fade2D_x64_${VS_VERSION}_Release.dll"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${CMAKE_CURRENT_SOURCE_DIR}/../fade2d/x64/libgmp-10.dll"
    "$<TARGET_FILE_DIR:${PROJECT_NAME}>/libgmp-10.dll"
)

target_compile_definitions(${PROJECT_NAME} PUBLIC DSYSTEM_FONT_DIR=\"${SYSTEM_FONT_DIR}\")