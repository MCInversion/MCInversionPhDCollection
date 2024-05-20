cmake_minimum_required(VERSION 3.14)
Project(poly2tri)

set(CMAKE_CONFIGURATION_TYPES Debug Release)

file(GLOB poly2tri_Src CONFIGURE_DEPENDS "*.h" "*.cpp")
add_library(poly2tri ${poly2tri_Src})