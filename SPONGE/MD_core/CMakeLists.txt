add_library(MD_core OBJECT MD_core.cu MD_core.cuh)
set_target_properties(MD_core PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
message("MD core found")