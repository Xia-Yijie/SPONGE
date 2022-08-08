add_library(TI_core OBJECT TI_core.cuh TI_core.cu)
set_target_properties(TI_core PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
message("TI core found.")