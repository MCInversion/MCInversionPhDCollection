cmake_minimum_required(VERSION 3.14)
project(pmp_tests)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(GOOGLE_TEST_ROOT
    external/googletest/googletest
    CACHE STRING "Google Test root")

include_directories("${PROJECT_SOURCE_DIR}/../../${GOOGLE_TEST_ROOT}"
                    "${PROJECT_SOURCE_DIR}/../../${GOOGLE_TEST_ROOT}/include")

set(GOOGLE_TEST_SOURCES
    "${PROJECT_SOURCE_DIR}/../../${GOOGLE_TEST_ROOT}/src/gtest-all.cc"
    "${PROJECT_SOURCE_DIR}/../../${GOOGLE_TEST_ROOT}/src/gtest_main.cc")

add_library(googletest ${GOOGLE_TEST_SOURCES})
file(GLOB UnitTests_Src "*.h" "*.cpp")

# Now simply link against gtest or gtest_main as needed. Eg
add_executable(${PROJECT_NAME} ${UnitTests_Src})
add_dependencies(${PROJECT_NAME} googletest)

target_include_directories(${PROJECT_NAME} PUBLIC
	"${CMAKE_SOURCE_DIR}/src/")

target_link_libraries(${PROJECT_NAME}
	googletest		 # gtest lib
	pmp              # tested lib
)

add_test(NAME pmp_TestSuite COMMAND ${PROJECT_NAME})

# ----- make sure CMAKE_CURRENT_SOURCE_DIR is readable in source ---------
#
# VS Test Explorer automatically sets working directory as CMAKE_CURRENT_BINARY_DIR,
# so we use DROOT_DIR in source files instead of std::filesystem::::current_path()
#
set(ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../..")
#message("ROOT_DIR: ${ROOT_DIR}")
target_compile_definitions(${PROJECT_NAME} PUBLIC DROOT_DIR=\"${ROOT_DIR}\")