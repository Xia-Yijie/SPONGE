# Flags
## Check CUDA Version, architecture and set compiler flags
include(FindCUDA/select_compute_arch)
CUDA_DETECT_INSTALLED_GPUS(INSTALLED_GPU_CCS_1)
string(STRIP "${INSTALLED_GPU_CCS_1}" INSTALLED_GPU_CCS_2)
string(REPLACE " " ";" INSTALLED_GPU_CCS_3 "${INSTALLED_GPU_CCS_2}")
string(REPLACE "." "" CUDA_ARCH_LIST "${INSTALLED_GPU_CCS_3}")
message("Supporting CUDA archs: \n ${CUDA_ARCH_LIST}")
list(GET CUDA_ARCH_LIST -1 CUDA_ARCH_LATEST)
message("Using architecture ${CUDA_ARCH_LATEST}")
message("The version of CUDA compiler is ${CMAKE_CUDA_COMPILER_VERSION}.")
if (CMAKE_CUDA_COMPILER_VERSION VERSION_LESS 11.0)
    message("Compiler version is less than 11.0, using flags:")
    message("\t -std=c++14 -O4 -rdc=true -lcudadevrt -lcufft -DNDEBUG --use_fast_math")
    add_compile_options(
        -lcudadevrt
        -lcufft
        -O4 
        --use_fast_math
        -DNDEBUG
    )
else()
    message("Compiler version is greater or equal to 11.0, using flags:")
    message("\t -std=c++14 -O4 -rdc=true -lcudadevrt -lcufft -DNDEBUG --use_fast_math")
    add_compile_options(
        -lcudadevrt
        -lcufft
        -O4 
        --use_fast_math
        -DNDEBUG
    )
endif()
# set cuda architectue
set(CMAKE_CUDA_ARCHITECTURES ${CUDA_ARCH_LATEST})

# Library
file(GLOB_RECURSE CURRENT_HEADERS CMAKE_CONFIGURE_DEPENDS *.cuh)
file(GLOB_RECURSE CURRENT_SOURCES CMAKE_CONFIGURE_DEPENDS *.cu)
list(FILTER CURRENT_HEADERS INCLUDE REGEX "${CMAKE_CURRENT_SOURCE_DIR}/.*/.*.cuh")
list(FILTER CURRENT_SOURCES INCLUDE REGEX "${CMAKE_CURRENT_SOURCE_DIR}/.*/.*.cu")
add_library(sponge STATIC ${CURRENT_HEADERS} ${CURRENT_SOURCES}) # use STATIC library for better performance
set_target_properties(sponge PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

# Executables
add_executable(SPONGE main.cuh main.cu common.cu common.cuh control.cu control.cuh)
set_target_properties(SPONGE PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
add_executable(SPONGE_NOPBC main_nopbc.cuh main_nopbc.cu common.cu common.cuh control.cu control.cuh)
set_target_properties(SPONGE_NOPBC PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
add_executable(SPONGE_FEP main_fep.cuh main_fep.cu common.cu common.cuh control.cu control.cuh)
set_target_properties(SPONGE_FEP PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
add_executable(SPONGE_TI main_ti.cuh main_ti.cu common.cu common.cuh control.cu control.cuh)
set_target_properties(SPONGE_TI PROPERTIES CUDA_SEPARABLE_COMPILATION ON)


# Link
target_link_libraries(SPONGE PUBLIC sponge -lcufft)
target_link_libraries(SPONGE_NOPBC PUBLIC sponge -lcufft)
target_link_libraries(SPONGE_FEP PUBLIC sponge -lcufft)
target_link_libraries(SPONGE_TI PUBLIC sponge -lcufft)