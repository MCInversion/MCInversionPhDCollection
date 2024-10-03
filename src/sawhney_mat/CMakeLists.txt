cmake_minimum_required(VERSION 3.14)
Project(sawhney_mat)

set(CMAKE_CONFIGURATION_TYPES Debug Release)

file(GLOB sawhney_mat_Src CONFIGURE_DEPENDS "*.h" "*.cpp")
add_library(${PROJECT_NAME} ${sawhney_mat_Src})