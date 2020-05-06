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
add_library(panel1x2_test STATIC panel1x2_test.f90)
target_link_libraries(panel1x2_test classdef naturalfruit )

add_library(panel1x3_test STATIC panel1x3_test.f90)
target_link_libraries(panel1x3_test classdef naturalfruit )

add_library(panel1x3NegPitch_test STATIC panel1x3NegPitch_test.f90)
target_link_libraries(panel1x3NegPitch_test classdef naturalfruit )

# Executables
add_executable(panel1x2_testrunner panel1x2_testrunner.f90)
target_link_libraries(panel1x2_testrunner panel1x2_test)

add_executable(panel1x3_testrunner panel1x3_testrunner.f90)
target_link_libraries(panel1x3_testrunner panel1x3_test)

add_executable(panel1x3NegPitch_testrunner panel1x3NegPitch_testrunner.f90)
target_link_libraries(panel1x3NegPitch_testrunner panel1x3NegPitch_test)

# Testing
add_test(NAME panel1x2 COMMAND panel1x2_testrunner)
add_test(NAME panel1x3 COMMAND panel1x3_testrunner)
add_test(NAME panel1x3NegPitch COMMAND panel1x3NegPitch_testrunner)
