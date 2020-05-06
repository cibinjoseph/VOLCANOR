cmake_minimum_required (VERSION 3.5)
project(VOLCANOR_TESTS Fortran)

# Two ways to include modules to the search path
# 1. Does not copy module files to bin/
include_directories(${CMAKE_BINARY_DIR})
# 2. Copies all module files to bin/. Clean up is difficult.
# set(CMAKE_Fortran_MODULE_DIRECTORY  ${CMAKE_BINARY_DIR})

# Common library
add_library(naturalfruit STATIC naturalfruit.f90)

# Libraries
add_library(panel2_test STATIC panel2_test.f90)
target_link_libraries(panel2_test classdef naturalfruit )

add_library(panel3_test STATIC panel3_test.f90)
target_link_libraries(panel3_test classdef naturalfruit )

# Executables
add_executable(panel2_testrunner panel2_testrunner.f90)
target_link_libraries(panel2_testrunner panel2_test)

add_executable(panel3_testrunner panel3_testrunner.f90)
target_link_libraries(panel3_testrunner panel3_test)

# Testing
add_test(NAME panel2 COMMAND panel2_testrunner)
add_test(NAME panel3 COMMAND panel3_testrunner)
