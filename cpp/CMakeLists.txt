cmake_minimum_required(VERSION 3.0)
project(LBM-ArrayFire
  VERSION 0.1.2
  LANGUAGES CXX)

find_package(ArrayFire)

if(ArrayFire_CPU_FOUND)
  # D2Q9 BGK
  add_executable(d2q9_bgk_channel_cpu d2q9_bgk_channel.cpp)
  target_link_libraries(d2q9_bgk_channel_cpu ArrayFire::afcpu)

  add_executable(d2q9_bgk_channel_packed_cpu d2q9_bgk_channel_packed.cpp)
  target_link_libraries(d2q9_bgk_channel_packed_cpu ArrayFire::afcpu)

  add_executable(d2q9_bgk_lid_cpu d2q9_bgk_lid.cpp)
  target_link_libraries(d2q9_bgk_lid_cpu ArrayFire::afcpu)

  # D3Q27 BGK
  add_executable(d3q27_bgk_channel_cpu d3q27_bgk_channel.cpp)
  target_link_libraries(d3q27_bgk_channel_cpu ArrayFire::afcpu)

  # add_executable(d3q27_bgk_lid_cpu d3q27_bgk_lid.cpp)
  # target_link_libraries(d3q27_bgk_lid_cpu ArrayFire::afcpu)

  # D3Q27 MRT
  add_executable(lbm_mrt_cpu lbm_mrt.cpp)
  target_link_libraries(lbm_mrt_cpu ArrayFire::afcpu)
endif()

if(ArrayFire_CUDA_FOUND)
  # D2Q9 BGK
  add_executable(d2q9_bgk_channel_cuda d2q9_bgk_channel.cpp)
  target_link_libraries(d2q9_bgk_channel_cuda ArrayFire::afcuda)

  add_executable(d2q9_bgk_channel_packed_cuda d2q9_bgk_channel_packed.cpp)
  target_link_libraries(d2q9_bgk_channel_packed_cuda ArrayFire::afcuda)

  add_executable(d2q9_bgk_lid_cuda d2q9_bgk_lid.cpp)
  target_link_libraries(d2q9_bgk_lid_cuda ArrayFire::afcuda)

  # D3Q27 BGK
  add_executable(d3q27_bgk_channel_cuda d3q27_bgk_channel.cpp)
  target_link_libraries(d3q27_bgk_channel_cuda ArrayFire::afcuda)

  # add_executable(d3q27_bgk_lid_cuda d3q27_bgk_lid.cpp)
  # target_link_libraries(d3q27_bgk_lid_cuda ArrayFire::afcuda)

  # D3Q27 MRT
  add_executable(lbm_mrt_cuda lbm_mrt.cpp)
  target_link_libraries(lbm_mrt_cuda ArrayFire::afcuda)
endif()

if(ArrayFire_OpenCL_FOUND)
  # D2Q9 BGK
  add_executable(d2q9_bgk_channel_opencl d2q9_bgk_channel.cpp)
  target_link_libraries(d2q9_bgk_channel_opencl ArrayFire::afopencl)

  add_executable(d2q9_bgk_channel_packed_opencl d2q9_bgk_channel_packed.cpp)
  target_link_libraries(d2q9_bgk_channel_packed_opencl ArrayFire::afopencl)

  add_executable(d2q9_bgk_lid_opencl d2q9_bgk_lid.cpp)
  target_link_libraries(d2q9_bgk_lid_opencl ArrayFire::afopencl)

  # D3Q27 BGK
  add_executable(d3q27_bgk_channel_opencl d3q27_bgk_channel.cpp)
  target_link_libraries(d3q27_bgk_channel_opencl ArrayFire::afopencl)

  # add_executable(d3q27_bgk_lid_opencl d3q27_bgk_lid.cpp)
  # target_link_libraries(d3q27_bgk_lid_opencl ArrayFire::afopencl)

  # D3Q27 MRT
  add_executable(lbm_mrt_opencl lbm_mrt.cpp)
  target_link_libraries(lbm_mrt_opencl ArrayFire::afopencl)
endif()
