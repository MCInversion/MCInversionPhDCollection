cmake_minimum_required(VERSION 3.14)

project(STBUtils)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

file(GLOB STBUtils_Src CONFIGURE_DEPENDS "*.h" "*.cpp")

add_library(${PROJECT_NAME} ${STBUtils_Src})
