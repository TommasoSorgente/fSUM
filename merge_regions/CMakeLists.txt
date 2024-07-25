cmake_minimum_required(VERSION 3.5)

project(merge_regions LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(merge_regions main.cpp)

set(cinolib_DIR "${PROJECT_SOURCE_DIR}/../cinolib")

set(CINOLIB_USES_OPENGL_GLFW_IMGUI ON)
set(CINOLIB_USES_VTK ON)

find_package(cinolib REQUIRED)

find_package(Boost REQUIRED)
add_definitions(-DCINOLIB_USES_BOOST)

if (APPLE)
    target_link_libraries(${PROJECT_NAME} PUBLIC cinolib /usr/local/lib/libshp.a)
elseif (UNIX)
    target_link_libraries(${PROJECT_NAME} PUBLIC cinolib /usr/local/lib/libshp.so)
endif()

add_compile_definitions(HOME_PATH="${PROJECT_SOURCE_DIR}/")

include(GNUInstallDirs)
install(TARGETS merge_regions
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
