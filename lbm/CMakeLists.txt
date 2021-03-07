cmake_minimum_required(VERSION 3.0)
project(LBM-ArrayFire
  VERSION 3.5.0
  LANGUAGES CXX)

find_package(ArrayFire)

# LBM BGK
if(ArrayFire_CPU_FOUND)
  add_executable(lbm_bgk_cpu lbm_bgk.cpp)
  target_link_libraries(lbm_bgk_cpu ArrayFire::afcpu)
endif()

if(ArrayFire_CUDA_FOUND)
  add_executable(lbm_bgk_cuda lbm_bgk.cpp)
  target_link_libraries(lbm_bgk_cuda ArrayFire::afcuda)
endif()

if(ArrayFire_OpenCL_FOUND)
  add_executable(lbm_bgk_opencl lbm_bgk.cpp)
  target_link_libraries(lbm_bgk_opencl ArrayFire::afopencl)
endif()

# LBM MRT
if(ArrayFire_CPU_FOUND)
  add_executable(lbm_mrt_cpu lbm_mrt.cpp)
  target_link_libraries(lbm_mrt_cpu ArrayFire::afcpu)
endif()

if(ArrayFire_CUDA_FOUND)
  add_executable(lbm_mrt_cuda lbm_mrt.cpp)
  target_link_libraries(lbm_mrt_cuda ArrayFire::afcuda)
endif()

if(ArrayFire_OpenCL_FOUND)
  add_executable(lbm_mrt_opencl lbm_mrt.cpp)
  target_link_libraries(lbm_mrt_opencl ArrayFire::afopencl)
endif()
