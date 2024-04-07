cmake_minimum_required(VERSION 3.14)
Project(Geometry)

set(CMAKE_CONFIGURATION_TYPES Debug Release)

file(GLOB Geometry_Src CONFIGURE_DEPENDS "*.h" "*.cpp")
add_library(Geometry ${Geometry_Src})

target_link_libraries(${PROJECT_NAME} pmp Utils quickhull)