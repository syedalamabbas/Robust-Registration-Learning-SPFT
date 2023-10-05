# Copyright (c) 2018, ArrayFire
# All rights reserved.
#
# This file is distributed under 3-clause BSD license.
# The complete license agreement can be obtained at:
# http://arrayfire.com/licenses/BSD-3-Clause

cmake_minimum_required(VERSION 3.0)
project(ArrayFire-Example-SPFT
  VERSION 3.5.0
  LANGUAGES CXX)

find_package(ArrayFire)

if(ArrayFire_CUDA_FOUND)
  add_executable(spft spft.cpp)
  target_link_libraries(spft ArrayFire::afcuda)
endif()
