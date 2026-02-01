cmake_minimum_required (VERSION 3.12)
project(VOLCANOR Fortran)

# Use variable FC for setting compiler before calling cmake
# FC=ifort cmake ..

# Default values for CMAKE_BUILD_TYPE (which will be empty initially)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE RELEASE)
  # message(WARNING "Setting CMAKE_BUILD_TYPE to: ${CMAKE_BUILD_TYPE}")
endif()

option(USE_BLAS_LAPACK "Link against LAPACK/BLAS" ON)
option(USE_OPENMP "Enable OpenMP parallelization" ON)

# Variables for gfortran
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  # Base common flags
  set(common_flags "-ffree-form -std=f2008 -fimplicit-none -ffree-line-length-none -fall-intrinsics")
  
  # Add LAPACK/BLAS if requested
  if(USE_BLAS_LAPACK)
    find_package(LAPACK)
    find_package(BLAS)
    if(NOT LAPACK_FOUND)
      message(WARNING "LAPACK might not be available")
    endif()
    if(NOT BLAS_FOUND)
      message(WARNING "BLAS might not be available")
    endif()
  endif()

  # Debug flags
  set(dbg_flags "-fbacktrace -O0 -Wall -Wextra -Wimplicit-interface -Wunused-parameter -Wcharacter-truncation -Wsurprising -Waliasing -fcheck=all -g -ffpe-trap=invalid,zero,overflow -fbounds-check")
  
  # Release flags
  set(rls_flags "-O2 -finline-small-functions -march=native -flto")
  if(USE_OPENMP)
    set(rls_flags "${rls_flags} -fopenmp")
  endif()
  
  set(COMPILERMSG "Compiler: GFORTRAN ${CMAKE_Fortran_COMPILER_VERSION}\n")
endif()

# Variables for ifort
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  # Base common flags
  set(common_flags "-free -implicitnone")
  
  # Add MKL if requested
  if(USE_BLAS_LAPACK)
    if (NOT DEFINED ENV{MKLROOT})
      message(WARNING "MKL might not be available")
    endif()
    set(common_flags "${common_flags} -qmkl")
  endif()
  
  # Set PROF default
  if(NOT DEFINED PROF)
    set(PROF FALSE)
  endif()
  
  # Debug flags
  if(PROF)
    set(dbg_flags "-traceback -O2 -g -pg -fpe0")
  else()
    set(dbg_flags "-traceback -O0 -warn all -check bounds -g -heap-arrays 200 -fpe0 -debug extended")
  endif()
  
  # Release flags
  if(PROF)
    set(rls_flags "-O3 -g -pg -fpe0 -xcore-avx2")
  else()
    set(rls_flags "-O3 -heap-arrays 200 -qopt-dynamic-align")
  endif()
  
  # Add OpenMP flags if requested
  if(USE_OPENMP)
    set(rls_flags "${rls_flags} -qopenmp -parallel")
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
if(USE_BLAS_LAPACK)
  add_library(Math STATIC ${PROJECT_SOURCE_DIR}/src/libMath.f90)
  target_link_libraries(Math PUBLIC ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
else()
  add_library(Math STATIC ${PROJECT_SOURCE_DIR}/src/libMath_noBLAS.f90)
endif()

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

# Target debug without LAPACK/BLAS
add_custom_target(volcanor_dbg_noBLAS
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug -DPROF=FALSE -DUSE_BLAS_LAPACK=OFF ${PROJECT_SOURCE_DIR}
  COMMENT ${COMPILERMSG}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target volcanor)

# Target release without LAPACK/BLAS
add_custom_target(volcanor_noBLAS
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release -DPROF=FALSE -DUSE_BLAS_LAPACK=OFF ${PROJECT_SOURCE_DIR}
  COMMENT ${COMPILERMSG}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target volcanor)

# Target release without LAPACK/BLAS and without OpenMP
add_custom_target(volcanor_noBLAS_noOMP
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release -DPROF=FALSE -DUSE_BLAS_LAPACK=OFF -DUSE_OPENMP=OFF ${PROJECT_SOURCE_DIR}
  COMMENT ${COMPILERMSG}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target volcanor)

# Target debug without LAPACK/BLAS and without OpenMP
add_custom_target(volcanor_dbg_noBLAS_noOMP
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug -DPROF=FALSE -DUSE_BLAS_LAPACK=OFF -DUSE_OPENMP=OFF ${PROJECT_SOURCE_DIR}
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
