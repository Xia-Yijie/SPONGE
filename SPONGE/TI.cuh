#ifndef MAIN_CUH
#define MAIN_CUH


#include "common.cuh"
#include "control.cuh"
#include "TI_core/TI_core.cuh"
//#include "kinetic_energy/kinetic_energy.cuh"
#include "bond/bond.cuh"
#include "bond_soft/bond_soft.cuh"
#include "angle/angle.cuh"
#include "dihedral/dihedral.cuh"
#include "LJ_soft_core/LJ_soft_core.cuh"
#include "nb14/nb14.cuh"
#include "neighbor_list/neighbor_list.cuh"
#include "PME_force/PME_force.cuh"

void Main_Initial(int argc, char *argv[]);
void Main_Iteration();
void Main_Print();
void Main_Calculation();
void Main_Volume_Update();
void Main_Clear();

#endif