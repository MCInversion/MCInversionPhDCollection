cmake_minimum_required(VERSION 3.14)
project(core_tests)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

file(GLOB core_UnitTests_Src "*.h" "*.cpp")

# Now simply link against gtest or gtest_main as needed. Eg
add_executable(${PROJECT_NAME} ${core_UnitTests_Src})
add_dependencies(${PROJECT_NAME} googletest)

target_include_directories(${PROJECT_NAME} PUBLIC
	"${CMAKE_SOURCE_DIR}/src/")

target_link_libraries(${PROJECT_NAME}
	googletest     # gtest lib
	Core           # tested lib
)

add_test(NAME core_TestSuite COMMAND ${PROJECT_NAME})

target_compile_definitions(${PROJECT_NAME} PUBLIC DROOT_DIR=\"${ROOT_DIR}\")