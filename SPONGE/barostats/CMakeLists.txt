add_library(barostats OBJECT andersen_barostat.cu andersen_barostat.cuh Berendsen_barostat.cu Berendsen_barostat.cuh MC_barostat.cu MC_barostat.cuh)
set_target_properties(barostats PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
message("Barostats found.")