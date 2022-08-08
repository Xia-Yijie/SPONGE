add_library(nb14 OBJECT nb14.cu nb14.cuh)
set_target_properties(nb14 PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
message("NB14 found.")