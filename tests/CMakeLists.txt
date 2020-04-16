cmake_minimum_required (VERSION 3.5)
project(VOLCANOR_TESTS Fortran)

# Common library
add_library(naturalfruit STATIC naturalfruit.f90)

# Libraries
add_library(panel2_test STATIC panel2_test.f90)
target_link_libraries(panel2_test classdef naturalfruit )

# Executables
add_executable(panel2_testrunner panel2_testrunner.f90)
target_link_libraries(panel2_testrunner panel2_test)

# Testing
add_test(panel2 panel2_testrunner)

