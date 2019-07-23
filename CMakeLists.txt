cmake_minimum_required (VERSION 3.5)
project(VOLCANOR Fortran)

# Use variable FC for setting compiler before calling cmake
# FC=ifort cmake ..

# Variables for gfortran
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(common_flags "-ffree-form -std=f2008 -fimplicit-none -ffree-line-length-none")
  set(dbg_flags "-fbacktrace -O0 -Wall -Wextra -Wimplicit-interface -Wunused-parameter \
  -Wcharacter-truncation -Wsurprising -Waliasing -fcheck=all -g -ffpe-trap=invalid,zero,overflow,underflow -fbounds-check")
  set(rls_flags "-O2 -fopenmp")
endif()

# Variables for ifort
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(common_flags "-free -implicitnone -r8")
  set(dbg_flags "-traceback -O0 -warn all -check bounds -g -fpe0 -debug extended -heap-arrays 200")
  set(rls_flags "-fast -O3 -qopenmp -parallel -heap-arrays 200 \
  -ansi-alias -qopt-jump-tables='large' -xcore-avx2")
endif()

# Set compiler flags using variables
set(CMAKE_Fortran_FLAGS "${common_flags}")
set(CMAKE_Fortran_FLAGS_DEBUG "${dbg_flags}")
set(CMAKE_Fortran_FLAGS_RELEASE "${rls_flags}")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

# Libraries
add_library(Math STATIC ${CMAKE_SOURCE_DIR}/libMath.f90)

add_library(C81 STATIC ${CMAKE_SOURCE_DIR}/libC81.f90)
target_link_libraries(C81 Math)

add_library(classdef STATIC ${CMAKE_SOURCE_DIR}/classdef.f90)
target_link_libraries(classdef Math C81)

add_library(Common STATIC ${CMAKE_SOURCE_DIR}/libCommon.f90)
target_link_libraries(Common classdef C81 Math)

add_library(Postprocess STATIC ${CMAKE_SOURCE_DIR}/libPostprocess.f90)
target_link_libraries(Postprocess classdef)

# Executables
add_executable(main.out ${CMAKE_SOURCE_DIR}/main.f90)
target_link_libraries(main.out Common Postprocess)

add_executable(gridgen.out ${CMAKE_SOURCE_DIR}/gridgen.f90)
target_link_libraries(gridgen.out Common)

# Command for target init target 
add_custom_command(
  OUTPUT ${CMAKE_SOURCE_DIR}/Results
  COMMAND echo "Creating Results/"
  COMMAND mkdir -p ${CMAKE_SOURCE_DIR}/Results)

# Target init
add_custom_target(init ALL
  DEPENDS ${CMAKE_SOURCE_DIR}/Results
  )

# Target fileclean 
add_custom_target(fileclean
  COMMAND echo "Clearing Results/"
  COMMAND rm -f ${CMAKE_SOURCE_DIR}/Results/*
  COMMAND rm -f ${CMAKE_SOURCE_DIR}/status.txt)

# Target release/run to regenerate Makefile
add_custom_target(run
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMAND reset
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target fileclean
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target main.out
  )

# Target debug/run_dbg to regenerate Makefile
add_custom_target(run_dbg
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Build ${CMAKE_SOURCE_DIR}
  COMMAND reset
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target fileclean
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target main.out
  )
