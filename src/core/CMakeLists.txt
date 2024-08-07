cmake_minimum_required(VERSION 3.14)
Project(Core)

set(CMAKE_CONFIGURATION_TYPES Debug Release)

file(GLOB Core_Src CONFIGURE_DEPENDS "*.h" "*.cpp")
add_library(${PROJECT_NAME} ${Core_Src})

target_link_libraries(${PROJECT_NAME} 
	pmp 
	Geometry
	SDF
	Utils
)