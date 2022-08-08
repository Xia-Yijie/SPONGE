add_library(No_PBC OBJECT Coulomb_Force_No_PBC.cu Coulomb_Force_No_PBC.cuh generalized_Born.cu generalized_Born.cuh Lennard_Jones_force_No_PBC.cu Lennard_Jones_force_No_PBC.cuh)
set_target_properties(No_PBC PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
message("No_PBC found.")