cmake_minimum_required(VERSION 3.14)
Project(SDFGen)

set(CMAKE_CONFIGURATION_TYPES Debug Release)

file(GLOB SDFGen_Src CONFIGURE_DEPENDS "*.h" "*.cpp")
add_library(SDFGen ${SDFGen_Src})