cmake_minimum_required(VERSION 3.5)
project(FESA LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(cinolib_DIR "${PROJECT_SOURCE_DIR}/cinolib")

set(CINOLIB_USES_OPENGL_GLFW_IMGUI ON)
set(CINOLIB_USES_VTK ON)

set(cinolib_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/cinolib)
find_package(cinolib REQUIRED)

find_package(Boost)

find_package(SQLite3)

if(NOT Boost_FOUND AND MSVC)

    if(NOT DEFINED ${CMAKE_TOOLCHAIN_FILE})
        #if(NOT DEFINED ENV{VCPKG_ROOT})
            #if(WIN32)
                set(VCPKG_ROOT $ENV{HOMEDRIVE}$ENV{HOMEPATH}/vcpkg_cpptruths)
            #else()
            #    set(VCPKG_ROOT $ENV{HOME}/.vcpkg_cpptruths)
            #endif()
        #else()
        #    set(VCPKG_ROOT $ENV{VCPKG_ROOT})
        #endif()

        if(NOT EXISTS ${VCPKG_ROOT})
            message("Cloning vcpkg in ${VCPKG_ROOT}")
            execute_process(COMMAND git clone https://github.com/Microsoft/vcpkg.git ${VCPKG_ROOT})
            # If a reproducible build is desired (and potentially old libraries are # ok), uncomment the
            # following line and pin the vcpkg repository to a specific githash.
            # execute_process(COMMAND git checkout 745a0aea597771a580d0b0f4886ea1e3a94dbca6 WORKING_DIRECTORY ${VCPKG_ROOT})
        else()
            # The following command has no effect if the vcpkg repository is in a detached head state.
            message("Auto-updating vcpkg in ${VCPKG_ROOT}")
            execute_process(COMMAND git pull WORKING_DIRECTORY ${VCPKG_ROOT})
        endif()

        if(NOT EXISTS ${VCPKG_ROOT}/README.md)
            message(FATAL_ERROR "***** FATAL ERROR: Could not clone vcpkg *****")
        endif()

        if(WIN32)
            set(BOOST_INCLUDEDIR ${VCPKG_ROOT}/installed/x86-windows/include)
            set(VCPKG_EXEC ${VCPKG_ROOT}/vcpkg.exe)
            set(VCPKG_BOOTSTRAP ${VCPKG_ROOT}/bootstrap-vcpkg.bat)
        else()
            set(VCPKG_EXEC ${VCPKG_ROOT}/vcpkg)
            set(VCPKG_BOOTSTRAP ${VCPKG_ROOT}/bootstrap-vcpkg.sh)
        endif()

        if(NOT EXISTS ${VCPKG_EXEC})
            message("Bootstrapping vcpkg in ${VCPKG_ROOT}")
            execute_process(COMMAND ${VCPKG_BOOTSTRAP} WORKING_DIRECTORY ${VCPKG_ROOT})
        endif()

        if(NOT EXISTS ${VCPKG_EXEC})
            message(FATAL_ERROR "***** FATAL ERROR: Could not bootstrap vcpkg *****")
        endif()

        set(CMAKE_TOOLCHAIN_FILE ${VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake CACHE STRING "")

        #message(STATUS "***** Checking project third party dependencies in ${VCPKG_ROOT} *****")
        execute_process(
                    COMMAND ${VCPKG_EXEC} install libboost --triplet x64-windows --recurse
                    WORKING_DIRECTORY ${VCPKG_ROOT})

        set (Boost_INCLUDE_DIR "${VCPKG_ROOT}/packages/boost_x64-windows/include" )
        set (Boost_LIBRARY "${SQLite3_INCLUDE_DIR}/../lib/libboost.lib" )

    endif()
endif()

add_definitions(-DCINOLIB_USES_BOOST)

set (SHAPELIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/shapelib")
set (SHAPELIB_BUILD "${SHAPELIB_DIR}/build-${CMAKE_SYSTEM_NAME}")

file (MAKE_DIRECTORY ${SHAPELIB_BUILD})

if (APPLE)
    set (SHAPELIB_LIB "${SHAPELIB_BUILD}/libshp.a")
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

include_directories(${SHAPELIB_DIR})

add_executable(${PROJECT_NAME} main.cpp parameters.run ${SHAPELIB_LIB})

target_link_libraries(${PROJECT_NAME} PUBLIC cinolib ${SHAPELIB_LIB})

add_compile_definitions(HOME_PATH="${PROJECT_SOURCE_DIR}/")
