add_library(constrain OBJECT constrain.cu constrain.cuh SETTLE.cu SETTLE.cuh SHAKE.cu SHAKE.cuh simple_constrain.cu simple_constrain.cuh)
set_target_properties(constrain PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
message("Constrain found.")