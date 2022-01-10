#include "TI.cuh"

CONTROLLER controller;
TI_CORE TI_core;
BOND bondA;
BOND bondB;
BOND_SOFT bond_soft;
ANGLE angleA;
ANGLE angleB;
DIHEDRAL dihedralA;
DIHEDRAL dihedralB;
NON_BOND_14 nb14A;
NON_BOND_14 nb14B;
NEIGHBOR_LIST neighbor_list;
LJ_SOFT_CORE lj_soft;
Particle_Mesh_Ewald pme;
//KINETIC_ENERGY kinetic;

int main(int argc, char *argv[])
{
	Main_Initial(argc, argv);
	//getchar();
	for (TI_core.input.current_frame = 1; TI_core.input.current_frame <= TI_core.input.frame_numbers; ++TI_core.input.current_frame)
	{
		Main_Calculation();
		//getchar();
		Main_Print();
		//getchar();
		if (TI_core.input.current_frame != TI_core.input.frame_numbers)
		{
			Main_Iteration();
		}
	}
	TI_core.Print_dH_dlambda_Average_To_Screen_And_Result_File(controller.fcresult);
	Main_Clear();
	return 0;
}

void Main_Initial(int argc, char *argv[])
{
	controller.Initial(argc, argv);
	TI_core.Initial(&controller);

	neighbor_list.Initial(&controller, TI_core.atom_numbers, TI_core.box_length, TI_core.nb.cutoff, TI_core.nb.skin, TI_core.need_nbl);
	neighbor_list.Neighbor_List_Update(TI_core.crd, TI_core.nb.d_excluded_list_start, TI_core.nb.d_excluded_list, TI_core.nb.d_excluded_numbers, 1);
	
	nb14A.Initial(&controller, TI_core.nb14_pertubated, "nb14A");
	nb14B.Initial(&controller, TI_core.nb14_pertubated, "nb14B");
	bondA.Initial(&controller, TI_core.bond_pertubated, "bondA");
	bondB.Initial(&controller, TI_core.bond_pertubated, "bondB");
	bond_soft.Initial(&controller, TI_core.bond_pertubated);
	angleA.Initial(&controller, TI_core.angle_pertubated, "angleA");
	angleB.Initial(&controller, TI_core.angle_pertubated, "angleB");
	dihedralA.Initial(&controller, TI_core.dihedral_pertubated, "dihedralA");
	dihedralB.Initial(&controller, TI_core.dihedral_pertubated, "dihedralB");

	//kinetic.Initial(&controller, TI_core.mass_pertubated);	
	lj_soft.Initial(&controller, TI_core.nb.cutoff, TI_core.box_length, TI_core.d_subsys_division,TI_core.lj_pertubated);
	pme.Initial(&controller, TI_core.atom_numbers, TI_core.box_length, TI_core.nb.cutoff, TI_core.d_charge_B_A, TI_core.charge_pertubated);

	controller.Input_Check();
	controller.Print_First_Line_To_Fcout();
}

void Main_Calculation()
{
	TI_core.data.dH_dlambda_current_frame = 0.0;
	TI_core.data.bondA_ene = bondA.Get_Energy(TI_core.uint_crd, TI_core.pbc.uint_dr_to_dr_cof);
	TI_core.data.bondB_ene = bondB.Get_Energy(TI_core.uint_crd, TI_core.pbc.uint_dr_to_dr_cof);
	TI_core.data.angleA_ene = angleA.Get_Energy(TI_core.uint_crd, TI_core.pbc.uint_dr_to_dr_cof);
	TI_core.data.angleB_ene = angleB.Get_Energy(TI_core.uint_crd, TI_core.pbc.uint_dr_to_dr_cof);
	TI_core.data.dihedralA_ene = dihedralA.Get_Energy(TI_core.uint_crd, TI_core.pbc.uint_dr_to_dr_cof);
	TI_core.data.dihedralB_ene = dihedralB.Get_Energy(TI_core.uint_crd, TI_core.pbc.uint_dr_to_dr_cof);
	TI_core.data.nb14A_EE_ene = nb14A.Get_14_CF_Energy(TI_core.uint_crd, TI_core.pbc.uint_dr_to_dr_cof);
	TI_core.data.nb14A_LJ_ene = nb14A.Get_14_LJ_Energy(TI_core.uint_crd, TI_core.pbc.uint_dr_to_dr_cof);
	TI_core.data.nb14B_EE_ene = nb14B.Get_14_CF_Energy(TI_core.uint_crd, TI_core.pbc.uint_dr_to_dr_cof);
	TI_core.data.nb14B_LJ_ene = nb14B.Get_14_LJ_Energy(TI_core.uint_crd, TI_core.pbc.uint_dr_to_dr_cof);

	TI_core.data.bond_soft_dH_dlambda = bond_soft.Get_Partial_H_Partial_Lambda(TI_core.uint_crd, TI_core.pbc.uint_dr_to_dr_cof);
	TI_core.data.lj_soft_dH_dlambda = lj_soft.Get_Partial_H_Partial_Lambda_With_Columb_Direct(TI_core.uint_crd, TI_core.d_charge, neighbor_list.d_nl, TI_core.d_charge_B_A, TI_core.charge_pertubated);
	if (lj_soft.is_initialized)
		TI_core.data.coul_direct_dH_dlambda = *lj_soft.h_sigma_of_dH_dlambda_direct;
	else
		TI_core.data.coul_direct_dH_dlambda = 0.0f;
	TI_core.data.lj_soft_long_range_correction = lj_soft.Partial_H_Partial_Lambda_Long_Range_Correction();
	TI_core.data.pme_dH_dlambda = pme.Get_Partial_H_Partial_Lambda_Without_Direct(TI_core.uint_crd, neighbor_list.d_nl, TI_core.d_charge, TI_core.pbc.uint_dr_to_dr_cof,
		TI_core.nb.d_excluded_list_start, TI_core.nb.d_excluded_list, TI_core.nb.d_excluded_numbers, TI_core.lj_pertubated);

	if (pme.is_initialized)
	{
		TI_core.data.pme_self_dH_dlambda = pme.cross_self_ene;
		TI_core.data.pme_reci_dH_dlambda = pme.cross_reciprocal_ene;
		TI_core.data.pme_corr_dH_dlambda = pme.cross_correction_ene;
		TI_core.data.coul_direct_dH_dlambda += pme.cross_direct_ene;
	}
	else
	{
		TI_core.data.pme_self_dH_dlambda = 0.0f;
		TI_core.data.pme_reci_dH_dlambda = 0.0f;
		TI_core.data.pme_corr_dH_dlambda = 0.0f;
	}
	//TI_core.data.kinetic_dH_dlambda = kinetic.Get_Partial_H_Partial_Lambda(TI_core.vel);
	
	TI_core.data.Sum_One_Frame();
}

void Main_Print()
{
	controller.Step_Print("frame", TI_core.input.current_frame);
	controller.Step_Print("dH_dlambda", TI_core.data.dH_dlambda_current_frame);
	controller.Step_Print("bondA", TI_core.data.bondA_ene);
	controller.Step_Print("bondB", TI_core.data.bondB_ene);
	controller.Step_Print("angleA", TI_core.data.angleA_ene);
	controller.Step_Print("angleB", TI_core.data.angleB_ene);
	controller.Step_Print("dihedralA", TI_core.data.dihedralA_ene);
	controller.Step_Print("dihedralB", TI_core.data.dihedralB_ene);
	controller.Step_Print("nb14A_LJ", TI_core.data.nb14A_LJ_ene);
	controller.Step_Print("nb14B_LJ", TI_core.data.nb14B_LJ_ene);
	controller.Step_Print("nb14A_EE", TI_core.data.nb14A_EE_ene);
	controller.Step_Print("nb14B_EE", TI_core.data.nb14B_EE_ene);
	controller.Step_Print("bond_soft", TI_core.data.bond_soft_dH_dlambda);
	controller.Step_Print("LJ(sc.)", TI_core.data.lj_soft_dH_dlambda);
	controller.Step_Print("Coul(direct.)", TI_core.data.coul_direct_dH_dlambda);
	controller.Step_Print("LR_corr", TI_core.data.lj_soft_long_range_correction);
	//controller.Step_Print("PME", TI_core.data.pme_dH_dlambda);
	controller.Step_Print("PME(reci.)", TI_core.data.pme_reci_dH_dlambda);
	controller.Step_Print("PME(self.)", TI_core.data.pme_self_dH_dlambda);
	controller.Step_Print("PME(corr.)", TI_core.data.pme_corr_dH_dlambda);
	//controller.Step_Print("kinetic", TI_core.data.kinetic_dH_dlambda);

	controller.Print_To_Screen_And_Fcout();
}

void Main_Volume_Update()
{
	lj_soft.Update_Volume(TI_core.box_length);
	pme.Update_Volume(TI_core.volume_change_factor);
	neighbor_list.Update_Volume(TI_core.box_length);
}

void Main_Iteration()
{
	TI_core.Read_Next_Frame();
	Main_Volume_Update();
	neighbor_list.Neighbor_List_Update(TI_core.crd, TI_core.nb.d_excluded_list_start, TI_core.nb.d_excluded_list, TI_core.nb.d_excluded_numbers, 1);
}

void Main_Clear()
{
	bondA.Clear();
	bondB.Clear();
	bond_soft.Clear();
	angleA.Clear();
	angleB.Clear();
	dihedralA.Clear();
	dihedralB.Clear();
	nb14A.Clear();
	nb14B.Clear();
	neighbor_list.Clear();
	lj_soft.Clear();
	pme.Clear();
	//kinetic.Clear();
	TI_core.Clear();
	controller.Clear();
}
