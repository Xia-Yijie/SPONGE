# CMAKE version required
cmake_minimum_required(VERSION 3.18)

# Set project
project(SPONGE-Project 
    VERSION 1.2
    DESCRIPTION "Simulation Package tOward Next GEneration molecular modelling"
    LANGUAGES C CXX CUDA)
set(CMAKE_VERBOSE_MAKEFILE ON)
# Check CUDA
set(CMAKE_CUDA_STANDARD 11)
if (CMAKE_CUDA_COMPILER)
    add_subdirectory(SPONGE)
else()
    message("CUDA compiler not found. Configure terminated.")
endif()





