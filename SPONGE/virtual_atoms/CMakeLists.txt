add_library(virtual_atoms OBJECT virtual_atoms.cu virtual_atoms.cuh)
set_target_properties(virtual_atoms PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
message("Virtual atoms found.")