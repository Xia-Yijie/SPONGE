# Flags
#set(CUDA_SEPARABLE_COMPILATION ON)
add_compile_options(
    -arch=sm_50
    -std=c++11
    -rdc=true
    -lcudadevrt
    -lcufft
    -O4 
    --use_fast_math
    -DNDEBUG
)

# Libraries
set(MODULES "angle;barostats;bond;cmap;constrain;crd_molecular_map;dihedral;Lennard_Jones_force;nb14;neighbor_list;PME_force;restrain;thermostats;virtual_atoms;FEP_core;MD_core;No_PBC;TI_core")
mark_as_advanced(MODULES)
message("Finding modules.")
foreach(_mod ${MODULES})
    message(${_mod})
    add_subdirectory(${_mod})
endforeach()

# Executables
add_executable(SPONGE main.cuh main.cu common.cu common.cuh control.cu control.cuh)
set_target_properties(SPONGE PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
add_executable(SPONGE_NOPBC main_nopbc.cuh main_nopbc.cu common.cu common.cuh control.cu control.cuh)
set_target_properties(SPONGE_NOPBC PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
add_executable(SPONGE_FEP main_fep.cuh main_fep.cu common.cu common.cuh control.cu control.cuh)
set_target_properties(SPONGE_FEP PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
add_executable(SPONGE_TI main_ti.cuh main_ti.cu common.cu common.cuh control.cu control.cuh)
set_target_properties(SPONGE_TI PROPERTIES CUDA_SEPARABLE_COMPILATION ON)


# Link
target_link_libraries(SPONGE PUBLIC
    MD_core 
    bond 
    bond_soft
    angle
    Urey_Bradley_force
    dihedral
    improper_dihedral
    cmap
    nb14
    neighbor_list
    Lennard_Jones_force
    LJ_soft
    PME_force
    thermostats
    barostats
    restrain
    constrain
    virtual_atoms
    crd_molecular_map
    -lcufft
    )

target_link_libraries(SPONGE_NOPBC PUBLIC
    MD_core
    bond
    angle
    Urey_Bradley_force
    improper_dihedral
    dihedral
    cmap
    nb14
    thermostats
    restrain
    constrain
    virtual_atoms
    crd_molecular_map
    No_PBC
    -lcufft
)

target_link_libraries(SPONGE_FEP PUBLIC
    FEP_core
    bond
    bond_soft
    angle
    dihedral
    nb14
    neighbor_list
    LJ_soft
    PME_force
    -lcufft
)

target_link_libraries(SPONGE_TI PUBLIC
    TI_core
    bond
    bond_soft
    angle
    dihedral
    nb14
    neighbor_list
    LJ_soft
    PME_force
    -lcufft
)