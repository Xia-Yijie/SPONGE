add_library(restrain OBJECT restrain.cu restrain.cuh)
set_target_properties(restrain PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
message("Restrain found.")