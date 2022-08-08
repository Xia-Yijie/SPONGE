add_library(crd_molecular_map OBJECT crd_molecular_map.cu crd_molecular_map.cuh)
set_target_properties(crd_molecular_map PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
message("Molecular map found.")