add_library(PME_force OBJECT PME_force.cu PME_force.cuh)
set_target_properties(PME_force PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
message("PME force found.")