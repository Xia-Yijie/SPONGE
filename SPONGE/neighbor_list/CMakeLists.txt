add_library(neighbor_list OBJECT neighbor_list.cu neighbor_list.cuh)
set_target_properties(neighbor_list PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
message("Neighbor list found.")