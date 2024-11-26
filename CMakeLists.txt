cmake_minimum_required(VERSION 2.8.3)
project(stanley_controller)

set(CMAKE_BUILD_TYPE Debug)
message(STATUS "BUILD TYPE: " ${CMAKE_BUILD_TYPE})

add_compile_options(-std=c++14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -Wextra")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

set(CMAKE_CXX_CLANG_TIDY
  clang-tidy;
  -header-filter=${PROJECT_SOURCE_DIR};
  -warnings-as-errors=*;
)

find_package(glog REQUIRED)

# GTest
enable_testing()
include(GoogleTest)
find_package(GTest REQUIRED)

include_directories(${PROJECT_SOURCE_DIR})

add_library(cubic_spline_interpolation stanley_controller/common/cubic_spline_interpolation.cc)
target_link_libraries(cubic_spline_interpolation glog::glog)

add_executable(cubic_spline_interpolation_test stanley_controller/common/cubic_spline_interpolation_test.cc)
target_link_libraries(cubic_spline_interpolation_test cubic_spline_interpolation)

add_library(stanley_controller stanley_controller/stanley_controller.cc)
target_link_libraries(stanley_controller glog::glog)

add_executable(stanley_controller_test stanley_controller/stanley_controller_test.cc)
target_link_libraries(stanley_controller_test stanley_controller cubic_spline_interpolation)

add_executable(spline_test stanley_controller/common/spline_test.cc)
target_link_libraries(spline_test glog::glog)
