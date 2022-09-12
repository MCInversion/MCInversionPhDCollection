cmake_minimum_required(VERSION 3.14)
Project(SDF)

set(CMAKE_CONFIGURATION_TYPES Debug Release)

file(GLOB SDF_Src CONFIGURE_DEPENDS "*.h" "*.cpp")
add_library(SDF ${SDF_Src})

target_link_libraries(${PROJECT_NAME} pmp Geometry)