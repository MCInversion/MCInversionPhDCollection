cmake_policy(SET CMP0048 NEW)
project(MCInversionPhDCollection
  VERSION 1.0)
cmake_minimum_required(VERSION 3.16.3)
cmake_policy(SET CMP0072 NEW)

option(PMP_BUILD_VIS "Build the PMP visualization tools" ON)
option(PMP_INSTALL "Install the PMP library and headers" ON)

# set output paths
set(PROJECT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

if(PMP_BUILD_VIS)
  set(GLFW_SOURCE_DIR "external/glfw")
  set(GLEW_SOURCE_DIR "external/glew")
  set(IMGUI_SOURCE_DIR "external/imgui")
  set(STBI_SOURCE_DIR "external/stb_image")

  # setup STB Image (place *before* GLFW since GLFW has an old copy of
  # stb_image_write.h)
  include_directories(${STBI_SOURCE_DIR})
  add_subdirectory(${STBI_SOURCE_DIR})

  # Building only the GLFW lib
  set(BUILD_SHARED_LIBS
      OFF
      CACHE BOOL "")
  set(GLFW_BUILD_EXAMPLES
      OFF
      CACHE BOOL "")
  set(GLFW_BUILD_TESTS
      OFF
      CACHE BOOL "")
  set(GLFW_BUILD_DOCS
      OFF
      CACHE BOOL "")
  set(GLFW_INSTALL
      OFF
      CACHE BOOL "")

  # setup use of GLFW and GLEW
  if(NOT EMSCRIPTEN)
    add_subdirectory(${GLFW_SOURCE_DIR} ${GLEW_SOURCE_DIR})

    include_directories(${GLFW_SOURCE_DIR}/include ${GLFW_SOURCE_DIR}/deps
                        ${GLEW_SOURCE_DIR}/include)

    add_definitions(-DGLEW_STATIC)
    add_library(glew STATIC ${GLEW_SOURCE_DIR}/src/glew.c
                            ${GLEW_SOURCE_DIR}/include)
    set_property(TARGET glew PROPERTY POSITION_INDEPENDENT_CODE ON)
    target_link_libraries(glew ${GLFW_LIBRARIES})
  endif()

  # setup IMGUI
  include_directories(${IMGUI_SOURCE_DIR})
  add_subdirectory(${IMGUI_SOURCE_DIR})
endif(PMP_BUILD_VIS)

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
if(CMAKE_GENERATOR MATCHES "Visual Studio 17 2022")
    set(CMAKE_CXX_STANDARD 23)  # Use C++23 for Visual Studio 2022
    message(STATUS "Setting C++ Standard to 23 for Visual Studio 2022")
elseif(CMAKE_GENERATOR MATCHES "Visual Studio 16 2019")
    set(CMAKE_CXX_STANDARD 20)  # Use default setting for Visual Studio 2019
    message(STATUS "Setting C++ Standard to 20 for Visual Studio 2019")
else()
    set(CMAKE_CXX_STANDARD ${CXX_STANDARD_DEFAULT})  # Use default for other generators
endif()

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# set default compiler flags for both gcc and clang
set(COMMON_CXX_FLAGS "-pedantic -Wall -Wextra -Wshadow -Werror")

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
# which directories to process
# ....................................................
add_subdirectory(src/fade2d)
add_subdirectory(src/poly2tri)
add_subdirectory(src/pmp)
add_subdirectory(src/sdf)
add_subdirectory(src/geometry)
add_subdirectory(src/quickhull)
add_subdirectory(src/utils)
file(GLOB MCInversionPhDCollection_Src CONFIGURE_DEPENDS "src/*.h" "src/*.cpp")
add_executable(MCInversionPhDCollection ${MCInversionPhDCollection_Src})
target_link_libraries(MCInversionPhDCollection
    pmp
    SDF
    Geometry
)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set_target_properties(
    glew
    imgui
    pmp_vis
    rply
    stb_image
    pmp
PROPERTIES FOLDER "PMP_Lib")

# ----- make sure CMAKE_CURRENT_SOURCE_DIR is readable in source ---------
#
# VS Test Explorer automatically sets working directory as CMAKE_CURRENT_BINARY_DIR,
# so we use DROOT_DIR in source files instead of std::filesystem::::current_path()
#
set(ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(MCInversionPhDCollection PUBLIC DROOT_DIR=\"${ROOT_DIR}\")
# Set MCInversionPhDCollection the startup project
if (MSVC)
	set_property(DIRECTORY ${ROOT_DIR} PROPERTY VS_STARTUP_PROJECT MCInversionPhDCollection)
endif()
# =======================================================

# add clang-format target if program is found
include(clang-format)
