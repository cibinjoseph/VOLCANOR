cmake_minimum_required (VERSION 3.12)
project(VOLCANOR Fortran)

# Use variable FC for setting compiler before calling cmake
# FC=ifort cmake ..

# Default values for CMAKE_BUILD_TYPE (which will be empty initially)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE RELEASE)
  # message(WARNING "Setting CMAKE_BUILD_TYPE to: ${CMAKE_BUILD_TYPE}")
endif()

# Variables for gfortran
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  # Check for LAPACK and BLAS
  find_package(LAPACK)
  find_package(BLAS)
  if(NOT LAPACK_FOUND)
    message(WARNING "LAPACK might not be available")
  endif()
  if(NOT BLAS_FOUND)
    message(WARNING "BLAS might not be available")
  endif()

  set(common_flags "-ffree-form -std=f2008 -fimplicit-none -ffree-line-length-none -fall-intrinsics -llapack -lblas")
  set(dbg_flags "-fbacktrace -O0 -Wall -Wextra -Wimplicit-interface -Wunused-parameter -Wcharacter-truncation -Wsurprising -Waliasing -fcheck=all -g -ffpe-trap=invalid,zero,overflow -fbounds-check")
  set(rls_flags "-O2 -fopenmp")
  set(COMPILERMSG "Compiler: GFORTRAN ${CMAKE_Fortran_COMPILER_VERSION}\n")
endif()

# Variables for ifort
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  if (NOT DEFINED ENV{MKLROOT})
    message(WARNING "MKL might not be available")
  endif()
  set(common_flags "-free -implicitnone -qmkl")
  if(NOT DEFINED PROF)
    set(PROF FALSE)
  endif()
  if (PROF)
    set(dbg_flags "-traceback -O2 -g -pg -fpe0")
    # -qopt-report-phase=vec -qopt-report=2")
    set(rls_flags "-O3 -g -pg -fpe0 -qopenmp -xcore-avx2 -parallel")
    # -qopt-report-phase=vec -qopt-report=2")
  else()
    set(dbg_flags "-traceback -O0 -warn all -check bounds -g -heap-arrays 200 -fpe0 -debug extended")
    set(rls_flags "-fast -O3 -qopenmp -parallel -heap-arrays 200 -ansi-alias -qopt-jump-tables='default' -xcore-avx2 -ipo -qopt-subscript-in-range -unroll-aggressive -qopt-dynamic-align") 
    # The following flags were removed on shifting to 
    # Ubuntu 18.04 due to reduction in performance, however for large panel
    # numbers it seems necessary in order to avoid a segfault
    # -heap-arrays 200
  endif()
  set(COMPILERMSG "Compiler: IFORT ${CMAKE_Fortran_COMPILER_VERSION}
  Build type : ${CMAKE_BUILD_TYPE}
  Release flags : ${rls_flags}
  Debug flags   : ${dbg_flags}\n")
endif()

# Set compiler flags using variables
set(CMAKE_Fortran_FLAGS "${common_flags}")
set(CMAKE_Fortran_FLAGS_DEBUG "${dbg_flags}")
set(CMAKE_Fortran_FLAGS_RELEASE "${rls_flags}")
set(CMAKE_BUILD_PARALLEL_LEVEL 4)

# Set other variables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(LOG_FILE ${CMAKE_BINARY_DIR}/build.log)

# Write to log file
file(WRITE ${LOG_FILE} ${COMPILERMSG})

# Testing
enable_testing()
set(CMAKE_CTEST_ARGUMENTS -V --output-on-failure --stop-on-failure)
add_subdirectory(tests)

# Libraries
add_library(Math STATIC ${PROJECT_SOURCE_DIR}/src/libMath.f90)

add_library(C81 STATIC ${PROJECT_SOURCE_DIR}/src/libC81.f90)
target_link_libraries(C81 Math)

add_library(classdef STATIC ${PROJECT_SOURCE_DIR}/src/classdef.f90)
target_link_libraries(classdef Math C81)

add_library(Common STATIC ${PROJECT_SOURCE_DIR}/src/libCommon.f90)
target_link_libraries(Common classdef C81 Math)

add_library(Postprocess STATIC ${PROJECT_SOURCE_DIR}/src/libPostprocess.f90)
target_link_libraries(Postprocess classdef)

# Executables
add_executable(volcanor ${PROJECT_SOURCE_DIR}/src/main.f90)
target_link_libraries(volcanor Common Postprocess)

add_executable(gridgen ${PROJECT_SOURCE_DIR}/src/gridgen.f90)
target_link_libraries(gridgen Common)

# Target release/volcanor_parallel to regenerate Makefile
# add_custom_target(volcanor_parallel
#   COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release -DPROF=FALSE ${PROJECT_SOURCE_DIR}
#   COMMENT ${COMPILERMSG}
#   COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target volcanor)

# Target debug/volcanor_dbg for debugging
add_custom_target(volcanor_dbg
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug -DPROF=FALSE ${PROJECT_SOURCE_DIR}
  COMMENT ${COMPILERMSG}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target volcanor)

# Target profile/volcanor_prof for profiling
add_custom_target(volcanor_prof
  COMMAND ${CMAKE_COMMAND} -DPROF=TRUE ${PROJECT_SOURCE_DIR}
  COMMENT ${COMPILERMSG}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target volcanor)

# Target profile/volcanor_dbg_prof for debugging and profiling
add_custom_target(volcanor_dbg_prof
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug -DPROF=TRUE ${PROJECT_SOURCE_DIR}
  COMMENT ${COMPILERMSG}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target volcanor)

# Target to build all testrunners since make test does not seem to work
add_custom_target(testrunners
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target libMath_testrunner
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target wing1x2_testrunner
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target wing1x3_testrunner
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target wing1x3NegPitch_testrunner
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target rotor1x2_testrunner
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target rotor1x2Rev_testrunner
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target rotor1x2NegPitch_testrunner)
