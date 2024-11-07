cmake_minimum_required(VERSION 3.5)

project(misclassification LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CINOLIB_USES_OPENGL_GLFW_IMGUI    ON)
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

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME} PUBLIC cinolib)

add_compile_definitions(HOME_PATH="${PROJECT_SOURCE_DIR}/")
