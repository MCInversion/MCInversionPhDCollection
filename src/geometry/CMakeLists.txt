cmake_minimum_required(VERSION 3.14)
Project(Geometry)

set(CMAKE_CONFIGURATION_TYPES Debug Release)

#add_compile_options(/bigobj) # for some Eigen solvers

file(GLOB Geometry_Src CONFIGURE_DEPENDS "*.h" "*.cpp")
add_library(${PROJECT_NAME} ${Geometry_Src})

# Ensure that the Fade2D include directories are accessible
target_include_directories(${PROJECT_NAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../fade2d/include_fade2d")

# Link Geometry library with other libraries, including Fade2D
target_link_libraries(${PROJECT_NAME} pmp Utils quickhull fade2d poly2tri)

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