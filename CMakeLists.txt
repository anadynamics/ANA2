cmake_minimum_required(VERSION 3.5.0)
project(ANA2 VERSION 2.0.0 LANGUAGES CXX)
set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")
set(CMAKE_BUILD_TYPE Release) 
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

###################################################################
if(NOT WIN32)
  string(ASCII 27 Esc)
  set(BoldRed     "${Esc}[1;31m")
  set(BoldGreen   "${Esc}[1;32m")
  set(BoldYellow  "${Esc}[1;33m")
  set(BoldBlue    "${Esc}[1;34m")
  set(BoldMagenta "${Esc}[1;35m")
  set(BoldCyan    "${Esc}[1;36m")
  set(BoldWhite   "${Esc}[1;37m")
endif()
###################################################################

add_executable(ANA2 "")
target_sources(ANA2
    PRIVATE

        src/Cavity.cpp
        src/CGALUtils.cpp
        src/ConvexHull.cpp
        src/ConvexHullFunctions.cpp
        src/main.cpp
        src/MD.cpp
        src/Modes.cpp
        src/Molecule.cpp
        src/NDD.cpp
        src/NDDUtils.cpp
        src/PDB.cpp
        src/ProgramOptions.cpp
        src/PrimitivesUtils.cpp
        src/Read.cpp
        src/Static.cpp
        src/TextUtils.cpp
        src/Utils.cpp
        src/Write.cpp


        include/ANA/Cavity.hpp
        include/ANA/CGALUtils.hpp
        include/ANA/ConvexHull.hpp
        include/ANA/ConvexHullFunctions.hpp
        include/ANA/Includes.hpp
        include/ANA/MD.hpp
        include/ANA/Modes.hpp
        include/ANA/Molecule.hpp
        include/ANA/NDD.hpp
        include/ANA/NDDUtils.hpp
        include/ANA/Options.hpp
        include/ANA/PDB.hpp
        include/ANA/Primitives.hpp
        include/ANA/PrimitivesUtils.hpp
        include/ANA/ProgramOptions.hpp
        include/ANA/Read.hpp
        include/ANA/Static.hpp
        include/ANA/TextUtils.hpp
        include/ANA/Utils.hpp
        include/ANA/Write.hpp
)

message("${Magenta}Trying to automatically locate CGAL, gmp and mpfr libraries.${ColourReset}")
find_package(CGAL CONFIG REQUIRED)
target_link_libraries(ANA2 PRIVATE CGAL)
target_link_libraries(ANA2 PRIVATE gmp)
target_link_libraries(ANA2 PRIVATE mpfr)

message("${Magenta}Trying to automatically locate Boost.${ColourReset}")
find_package(Boost REQUIRED COMPONENTS program_options)
target_link_libraries(ANA2 PRIVATE boost_program_options)

message("${Magenta}Building own chemfiles library.${ColourReset}")
add_subdirectory(chemfiles EXCLUDE_FROM_ALL)
target_link_libraries(ANA2 PRIVATE chemfiles)

message("${Magenta}Using chemfiles's fmt header-only library.${ColourReset}")
include_directories(build/chemfiles/external/fmt)

target_include_directories(ANA2
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

#message("${Magenta}Including cotire for pre-compiled header.${ColourReset}")
#include(cotire)
#cotire(ANA2)

target_compile_options(ANA2 PRIVATE -Wall -Wextra -Wpedantic)
target_compile_features(ANA2 PRIVATE cxx_std_17)

option (FORCE_COLORED_OUTPUT "Always produce ANSI-colored output (GNU/Clang only)." ON)
if (${FORCE_COLORED_OUTPUT})
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        target_compile_options(ANA2 PRIVATE -fdiagnostics-color=always)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        target_compile_options(ANA2 PRIVATE -fcolor-diagnostics)
    endif ()
endif ()
