cmake_minimum_required (VERSION 3.12)
project(VOLCANOR_TESTS Fortran)

# Check for python3 for full unit tests
find_package(Python3 COMPONENTS Interpreter)

# Two ways to include modules to the search path
# 1. Does not copy module files to bin/
include_directories(${CMAKE_BINARY_DIR})
# 2. Copies all module files to bin/. Clean up is difficult.
# set(CMAKE_Fortran_MODULE_DIRECTORY  ${CMAKE_BINARY_DIR})

# Common library
add_library(naturalfruit STATIC naturalfruit.f90)

# Libraries
add_library(libMath_test STATIC libMath_test.f90)
target_link_libraries(libMath_test Math naturalfruit )

add_library(panel1x2_test STATIC panel1x2_test.f90)
target_link_libraries(panel1x2_test classdef naturalfruit )

add_library(panel1x3_test STATIC panel1x3_test.f90)
target_link_libraries(panel1x3_test classdef naturalfruit )

add_library(panel1x3NegPitch_test STATIC panel1x3NegPitch_test.f90)
target_link_libraries(panel1x3NegPitch_test classdef naturalfruit )

# Executables
add_executable(libMath_testrunner libMath_testrunner.f90)
target_link_libraries(libMath_testrunner libMath_test)

add_executable(panel1x2_testrunner panel1x2_testrunner.f90)
target_link_libraries(panel1x2_testrunner panel1x2_test)

add_executable(panel1x3_testrunner panel1x3_testrunner.f90)
target_link_libraries(panel1x3_testrunner panel1x3_test)

add_executable(panel1x3NegPitch_testrunner panel1x3NegPitch_testrunner.f90)
target_link_libraries(panel1x3NegPitch_testrunner panel1x3NegPitch_test)

# Unit testing
add_test(NAME libMath COMMAND libMath_testrunner)
add_test(NAME panel1x2 COMMAND panel1x2_testrunner)
add_test(NAME panel1x3 COMMAND panel1x3_testrunner)
add_test(NAME panel1x3NegPitch COMMAND panel1x3NegPitch_testrunner)

# Full-case testing
if (NOT DEFINED TESTFULL)
  set(TESTFULL TRUE)
endif()

if (TESTFULL)
  if (NOT Python3_Interpreter_FOUND)
    message(WARNING "Python 3 might not be available")
    message(WARNING "Full-case tests require Python 3")
    message(WARNING "Add -DTESTFULL=False to disable full-case tests")
  endif()
  add_test(NAME elevate 
    COMMAND python3 -m unittest discover ${CMAKE_CURRENT_SOURCE_DIR})
endif()
