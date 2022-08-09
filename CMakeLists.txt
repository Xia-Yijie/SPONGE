# CMAKE version required
cmake_minimum_required(VERSION 3.18)
set(CMAKE_VERBOSE_MAKEFILE ON)

# Check CUDA
include(CheckLanguage)
check_language(CUDA)
if (CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
else()
    message("Project will not be built.")
endif()

# Set project
project(SPONGE-Project 
    VERSION 1.2
    DESCRIPTION "Simulation Package tOward Next GEneration molecular modelling"
    LANGUAGES C CXX CUDA)

add_subdirectory(SPONGE)




