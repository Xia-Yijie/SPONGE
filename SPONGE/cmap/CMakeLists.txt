add_library(cmap OBJECT cmap.cu cmap.cuh)
set_target_properties(cmap PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
message("CMAP found.")