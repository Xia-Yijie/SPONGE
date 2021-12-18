#ifndef MAIN_CUH
#define MAIN_CUH


#include "common.cuh"
#include "control.cuh"
#include "MD_core/MD_core.cuh"
#include "bond/bond.cuh"
#include "angle/angle.cuh"
#include "dihedral/dihedral.cuh"
#include "nb14/nb14.cuh"
#include "neighbor_list/neighbor_list.cuh"
#include "Lennard_Jones_force/Lennard_Jones_force.cuh"
#include "PME_force/PME_force.cuh"
#include "Langevin_MD/Langevin_MD.cuh"
#include "Langevin_MD/LiuJian_MD.cuh"
#include "MC_barostat/MC_barostat.cuh"
#include "Berendsen_barostat/Berendsen_barostat.cuh"
#include "restrain/restrain.cuh"
#include "simple_constrain/simple_constrain.cuh"
#include "virtual_atoms/virtual_atoms.cuh"
#include "crd_molecular_map/crd_molecular_map.cuh"


void Main_Initial(int argc, char *argv[]);
void Main_Clear();

void Main_Calculate_Force();
void Main_Iteration();
void Main_Print();

void Main_Volume_Change(double factor);
void Main_Volume_Change_Largely();

#endif //MAIN_CUH