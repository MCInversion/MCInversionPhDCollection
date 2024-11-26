cmake_minimum_required(VERSION 3.14)
Project(Utils)

set(CMAKE_CONFIGURATION_TYPES Debug Release)

file(GLOB Utils_Src CONFIGURE_DEPENDS "*.h" "*.cpp" "nlohmann/*.hpp")
add_library(Utils ${Utils_Src})

# target_link_libraries(${PROJECT_NAME} pmp)