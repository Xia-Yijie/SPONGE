add_library(FEP_core OBJECT FEP_core.cu FEP_core.cuh)
set_target_properties(FEP_core PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
message("FEP core found.")