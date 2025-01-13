cmake_minimum_required(VERSION 3.5)
project(segmentation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CINOLIB_USES_OPENGL_GLFW_IMGUI    ON )
set(CINOLIB_USES_TETGEN               OFF)
set(CINOLIB_USES_TRIANGLE             OFF)
set(CINOLIB_USES_SHEWCHUK_PREDICATES  OFF)
set(CINOLIB_USES_INDIRECT_PREDICATES  OFF)
set(CINOLIB_USES_GRAPH_CUT            OFF)
set(CINOLIB_USES_BOOST                OFF)
set(CINOLIB_USES_VTK                  OFF) # optional, if added VTK formats will be supported everywhere
set(CINOLIB_USES_SPECTRA              OFF)
set(CINOLIB_USES_CGAL_GMP_MPFR        OFF)

set(cinolib_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../external/cinolib)
find_package(cinolib REQUIRED)

set(tclap_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../external/tclap/include)

set (SHAPELIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../external/shapelib")
set (SHAPELIB_BUILD "${SHAPELIB_DIR}/build-${CMAKE_SYSTEM_NAME}")

file (MAKE_DIRECTORY ${SHAPELIB_BUILD})

if (APPLE)
    set (SHAPELIB_LIB "${SHAPELIB_BUILD}/libshp.dylib")
elseif (UNIX)
    set (SHAPELIB_LIB "${SHAPELIB_BUILD}/libshp.so")
elseif (MSVC)
    set (SHAPELIB_LIB "${SHAPELIB_BUILD}/dll/Debug/shp.lib")
endif()

if(NOT EXISTS "${SHAPELIB_LIB}")
    message ("Building shapelib... ${SHAPELIB_LIB}")
    add_custom_command(
        OUTPUT "${SHAPELIB_LIB}"
        COMMAND cmake .. && cmake --build .
        WORKING_DIRECTORY "${SHAPELIB_BUILD}" )
endif()

include_directories(${SHAPELIB_DIR} ${tclap_DIR})

add_executable(${PROJECT_NAME} main.cpp ${SHAPELIB_LIB})

target_link_libraries(${PROJECT_NAME} PUBLIC cinolib ${SHAPELIB_LIB})

add_compile_definitions(HOME_PATH="${PROJECT_SOURCE_DIR}/")
