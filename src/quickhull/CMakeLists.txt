cmake_minimum_required(VERSION 3.14)
Project(quickhull)

set(CMAKE_CONFIGURATION_TYPES Debug Release)

file(GLOB quickhull_Src CONFIGURE_DEPENDS "*.h" "*.cpp" "*.hpp" "Structs/*.cpp" "Structs/*.hpp")
add_library(quickhull ${quickhull_Src})