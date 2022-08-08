add_library(thermostats OBJECT Andersen_thermostat.cu Andersen_thermostat.cuh Berendsen_thermostat.cu Berendsen_thermostat.cuh Langevin_MD.cu Langevin_MD.cuh Middle_Langevin_MD.cu Middle_Langevin_MD.cuh nose_hoover_chain.cu nose_hoover_chain.cuh)
set_target_properties(thermostats PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
message("Thermostats found.")