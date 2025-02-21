cmake_policy(SET CMP0048 NEW)

set(MAIN_PROJECT_NAME "MCInversionPhDCollection")

project(${MAIN_PROJECT_NAME} VERSION 1.0)
cmake_minimum_required(VERSION 3.16.3)
cmake_policy(SET CMP0072 NEW)

option(PMP_INSTALL "Install the PMP library and headers" ON)
option(BUILD_UNIT_TESTS "Build the unit testing projects" ON)

# set output paths
set(PROJECT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# setup Eigen
set(EIGEN_SOURCE_DIR "external/eigen")
include_directories(${EIGEN_SOURCE_DIR})

# setup PLY
set(RPLY_SOURCE_DIR "external/rply")
include_directories(${RPLY_SOURCE_DIR})
add_subdirectory(${RPLY_SOURCE_DIR})

# setup Nanoflann
set(NANOFLANN_SRC_DIR "external/nanoflann")
include_directories(${NANOFLANN_SRC_DIR})

# setup VCG
set(VCG_SRC_DIR "external/vcglib")
include_directories(${VCG_SRC_DIR})

# setup Spectra
# add_subdirectory(external/spectra)

# setup Tetgen
set(TETGEN_SRC_DIR "external/tetgen")
include_directories(${TETGEN_SRC_DIR})

# setup Triange (Libgl)
set(TRIANGLE_SRC_DIR "external/triangle")
add_subdirectory(${TRIANGLE_SRC_DIR})
include_directories(${TRIANGLE_SRC_DIR})

include(AddFileDependencies)
include_directories(${PROJECT_SOURCE_DIR}/src/)

# setup for code coverage testing
if(CMAKE_BUILD_TYPE STREQUAL "Debug"
   AND ENABLE_COVERAGE
   AND NOT WIN32)
  enable_testing()

  # set compiler flags
  set(CMAKE_CXX_FLAGS "-g -O0 --coverage")

  # find required tools
  find_program(LCOV lcov REQUIRED)
  find_program(GENHTML genhtml REQUIRED)
  find_package(OpenMP REQUIRED)

if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

  # add coverage target
  add_custom_target(
    coverage
    # gather data
    COMMAND
      ${LCOV} --directory . --capture --exclude '*/external/*' --exclude
      '*/tests/*' --exclude '/usr/*' --exclude '/Applications/*' --exclude
      'v1/*' --output-file coverage.info
    # generate report
    COMMAND ${GENHTML} --demangle-cpp -o coverage coverage.info
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
endif()

# Check the generator and adjust the C++ standard if necessary
set(CMAKE_CXX_STANDARD 20)  # Use default setting for Visual Studio 2019
message(STATUS "Setting C++ Standard to 20 for Visual Studio 2019")

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set default compiler flags for both gcc and clang
set(COMMON_CXX_FLAGS "-pedantic -Wall -Wextra -Wshadow -Werror")
# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /Od /Zi") # Disable optimization and add debug information in Release mode

if((UNIX OR APPLE) AND CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} ${COMMON_CXX_FLAGS} -Wno-deprecated-copy")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "7.0.0")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-int-in-bool-context")
  endif()
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND NOT EMSCRIPTEN)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMMON_CXX_FLAGS}")
endif()

if(WIN32)
  set(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} -D_USE_MATH_DEFINES -DNOMINMAX -D_CRT_SECURE_NO_WARNINGS"
  )
endif()

if(MSVC)
  add_compile_options(/MP)     # parallel compilation
  add_compile_options(/bigobj) # big object file
  add_compile_options(/wd4267) # conversion from size_t
  add_compile_options(/wd4244) # conversion from double to float
  add_compile_options(/wd4305) # truncation from double to float
endif()

if(EMSCRIPTEN)
  add_compile_options(--no-heap-copy)
  add_link_options(
    "SHELL:-s WASM=1 -s USE_WEBGL2=1 -s USE_GLFW=3 -s ALLOW_MEMORY_GROWTH=1 -s MINIFY_HTML=0 -s DISABLE_EXCEPTION_CATCHING=0"
  )
  set(CMAKE_EXECUTABLE_SUFFIX ".html")
endif()

# make a release build by default
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type specified. Using 'Release' as default.")
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE
        STRING
        "The build type. Possible values: Debug Release RelWithDebInfo MinSizeRel."
        FORCE)
endif()

if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb")
endif()

# choose size of the scalar type
if(PMP_SCALAR_TYPE EQUAL 64)
  message(STATUS "Using 64-bit scalar type")
  add_definitions(-DPMP_SCALAR_TYPE_64)
endif()

# choose size of the index type
if(PMP_INDEX_TYPE EQUAL 64)
  message(STATUS "Using 64-bit index type")
  add_definitions(-DPMP_INDEX_TYPE_64)
endif()

# setup clang-tidy if program found
option(WITH_CLANG_TIDY "Run clang-tidy checks" OFF)
include(clang-tidy)

# ====================================================
#         useful paths
# ....................................................

# ----- make sure CMAKE_CURRENT_SOURCE_DIR is readable in source ---------
#
# VS Test Explorer automatically sets working directory as CMAKE_CURRENT_BINARY_DIR,
# so we use DROOT_DIR in source files instead of std::filesystem::::current_path()
#
set(ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

if (CMAKE_SYSTEM_NAME STREQUAL "Windows")
    set(SYSTEM_FONT_DIR "C:/Windows/Fonts")
elseif (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(SYSTEM_FONT_DIR "/Library/Fonts")
elseif (CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(SYSTEM_FONT_DIR "/usr/share/fonts")
else()
    message(FATAL_ERROR "Unsupported operating system")
endif()

# ====================================================
# which directories to process
# ....................................................
add_subdirectory(src/core)
add_subdirectory(src/fade2d)
add_subdirectory(src/poly2tri)
add_subdirectory(src/pmp)
add_subdirectory(src/sdf)
add_subdirectory(src/geometry)
add_subdirectory(src/sawhney_mat)
add_subdirectory(src/stb_utils)
add_subdirectory(src/utils)

file(GLOB Experiment_Src CONFIGURE_DEPENDS "src/experiments/*.h" "src/experiments/*.cpp")
file(GLOB Main_Src CONFIGURE_DEPENDS "src/*.h" "src/*.cpp")

source_group("Experiments" FILES ${Experiment_Src})
source_group("Main" FILES ${Main_Src})

add_executable(${PROJECT_NAME} ${Main_Src} ${Experiment_Src})

target_link_libraries(${PROJECT_NAME}
    pmp
    SDF
    Geometry
    Core
    STBUtils
)

target_include_directories(${PROJECT_NAME} PUBLIC "${CMAKE_SOURCE_DIR}/src/")

target_compile_definitions(${PROJECT_NAME} PUBLIC 
    DROOT_DIR=\"${ROOT_DIR}\"
    #DSYSTEM_FONT_DIR=\"${SYSTEM_FONT_DIR}\"
)


if (BUILD_UNIT_TESTS)
    set(GOOGLE_TEST_ROOT
        external/googletest/googletest
        CACHE STRING "Google Test root")

    include_directories("${PROJECT_SOURCE_DIR}/${GOOGLE_TEST_ROOT}"
                        "${PROJECT_SOURCE_DIR}/${GOOGLE_TEST_ROOT}/include")

    set(GOOGLE_TEST_SOURCES
        "${PROJECT_SOURCE_DIR}/${GOOGLE_TEST_ROOT}/src/gtest-all.cc"
        "${PROJECT_SOURCE_DIR}/${GOOGLE_TEST_ROOT}/src/gtest_main.cc")
    add_library(googletest ${GOOGLE_TEST_SOURCES})

    add_subdirectory(src/pmp_tests)
    add_subdirectory(src/geometry_tests)
    add_subdirectory(src/sdf_tests)
    add_subdirectory(src/core_tests)
endif()

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set_target_properties(
    pmp
PROPERTIES FOLDER "PMP_Lib")

# Set PROJECT_NAME the startup project
if (MSVC)
	set_property(DIRECTORY ${ROOT_DIR} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME})
endif()
# =======================================================

# add clang-format target if program is found
include(clang-format)
