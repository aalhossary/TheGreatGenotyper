cmake_minimum_required(VERSION 3.1)
project(TheGreatGenotyper)

set (CMAKE_CXX_STANDARD 17)




add_executable(PopMergeVCF PopMergeVCF.cpp)

# Regression test: the output for a line must not depend on earlier lines.
enable_testing()
add_test(NAME PopMergeVCF_line_independence
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_line_independence.sh $<TARGET_FILE:PopMergeVCF>)
