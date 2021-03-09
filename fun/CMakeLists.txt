cmake_minimum_required(VERSION 3.0)
project(LBM-ArrayFire
  VERSION 0.1.1
  LANGUAGES CXX)

find_package(ArrayFire)

if(ArrayFire_CPU_FOUND)
  add_executable(d2q9_bgk_channel_packed_cpu d2q9_bgk_channel_packed.cpp)
  target_link_libraries(d2q9_bgk_channel_packed_cpu ArrayFire::afcpu)
endif()

if(ArrayFire_CUDA_FOUND)
  add_executable(d2q9_bgk_channel_packed_cuda d2q9_bgk_channel_packed.cpp)
  target_link_libraries(d2q9_bgk_channel_packed_cuda ArrayFire::afcuda)
endif()

if(ArrayFire_OpenCL_FOUND)
  add_executable(d2q9_bgk_channel_packed_opencl d2q9_bgk_channel_packed.cpp)
  target_link_libraries(d2q9_bgk_channel_packed_opencl ArrayFire::afopencl)
endif()