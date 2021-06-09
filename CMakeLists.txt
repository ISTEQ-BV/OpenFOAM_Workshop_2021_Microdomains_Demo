cmake_minimum_required(VERSION 3.10.0 FATAL_ERROR)
project(openfoam_workshop_2011_microdomains_demo VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

add_compile_options(-std=c++11 -Wall -Wextra -Wno-unused-parameter -Wno-invalid-offsetof -g)

if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "CMAKE_BUILD_TYPE is not specified, setting to 'Release'")
    set(CMAKE_BUILD_TYPE Release)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    message(STATUS "Building release configuration")
    message(STATUS "CMAKE_CXX_COMPILER_ID: ${CMAKE_CXX_COMPILER_ID}")
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        add_compile_options(-O3 -march=native)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        add_compile_options(-O3 -march=native -Wno-unknown-pragmas)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        add_compile_options(-xHOST -no-prec-div)
    else()
        message(FATAL_ERROR "Unknown compiler type: ${CMAKE_CXX_COMPILER_ID}")
    endif()
else()
    message(STATUS "Building debug configuration")
    add_compile_options(-g -ggdb3 -O0)
endif()

set(FOAM_SRC    "~/OpenFOAM/OpenFOAM-v2012/src" CACHE PATH "Path to OpenFOAM source directory")
set(FOAM_LIBBIN "~/OpenFOAM/OpenFOAM-v2012/platforms/linux64GccDPInt32Opt/lib" CACHE PATH "Path to OpenFOAM lib directory")

message(STATUS "FOAM_SRC: ${FOAM_SRC}")
message(STATUS "FOAM_LIBBIN: ${FOAM_LIBBIN}")

include_directories( SYSTEM
    ${FOAM_SRC}/OpenFOAM/lnInclude
    ${FOAM_SRC}/finiteVolume/lnInclude
    ${FOAM_SRC}/OSspecific/POSIX/lnInclude
    ${FOAM_SRC}/meshTools/lnInclude
    )
link_directories(
    ${FOAM_LIBBIN}
    ${FOAM_LIBBIN}/sys-openmpi
    )
add_definitions(-Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -DNoRepository)

find_package(MPI COMPONENTS CXX REQUIRED)

add_subdirectory(src)

enable_testing()
