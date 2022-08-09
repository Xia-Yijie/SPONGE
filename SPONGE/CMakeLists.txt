# Flags
## Check CUDA Version and set compiler flags
message("The version of CUDA compiler is ${CMAKE_CUDA_COMPILER_VERSION}.")
if (CMAKE_CUDA_COMPILER_VERSION VERSION_LESS 11.0)
    message("Compiler version is less than 11.0, using flags:")
    message("\t -O4 -std=c++11 -arch=sm_50 -rdc=true -lcudadevrt -lcufft -DNDEBUG --use_fast_math")
    add_compile_options(
        -arch=sm_50
        -std=c++11
        -rdc=true
        -lcudadevrt
        -lcufft
        -O4 
        --use_fast_math
        -DNDEBUG
    )
else()
    message("Compiler version is greater or equal to 11.0, using flags:")
    message("\t -O4 -std=c++14 -arch=sm_75 -rdc=true -lcudadevrt -lcufft -DNDEBUG --use_fast_math")
    add_compile_options(
        -arch=sm_75
        -std=c++14
        -rdc=true
        -lcudadevrt
        -lcufft
        -O4 
        --use_fast_math
        -DNDEBUG
    )
endif()

# Library
file(GLOB_RECURSE CURRENT_HEADERS CMAKE_CONFIGURE_DEPENDS *.cuh)
file(GLOB_RECURSE CURRENT_SOURCES CMAKE_CONFIGURE_DEPENDS *.cu)
list(FILTER CURRENT_HEADERS INCLUDE REGEX "${CMAKE_CURRENT_SOURCE_DIR}/.*/.*.cuh")
list(FILTER CURRENT_SOURCES INCLUDE REGEX "${CMAKE_CURRENT_SOURCE_DIR}/.*/.*.cu")
add_library(libsponge STATIC ${CURRENT_HEADERS} ${CURRENT_SOURCES}) # use STATIC library for better performance
set_target_properties(libsponge PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

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
target_link_libraries(SPONGE PUBLIC libsponge -lcufft)
target_link_libraries(SPONGE_NOPBC PUBLIC libsponge -lcufft)
target_link_libraries(SPONGE_FEP PUBLIC libsponge -lcufft)
target_link_libraries(SPONGE_TI PUBLIC libsponge -lcufft)