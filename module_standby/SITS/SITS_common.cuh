//#ifndef SITS_COMMON_CUH
//#define SITS_COMMON_CUH
//
//#include "..\common.cuh"
//#include "..\AMBER_file_IO\AMBER_file_IO.cuh"
//#include "..\neighbor_list\neighbor_list.cuh"
//#include "..\Lennard_Jones_force\Lennard_Jones_force.cuh"
//#include "..\MD_core\MD_core.cuh"
//#include "..\Langvein_MD\Langvein_MD.cuh"
//#include "..\Langvein_MD\LiuJian_MD.cuh"
//#include "..\bond_force\bond_force.cuh"
//#include "..\angle_force\angle_force.cuh"
//#include "..\dihedral_force\dihedral_force.cuh"
//#include "..\PME_force\PME_force.cuh"
//#include "float.h"
//
////本头文件是给SITS方法所使用的主头文件，依赖于普通MD的类容。
////使用方法见下面两篇文献
////A selective integrated tempering method
////Self-adaptive enhanced sampling in the energy and trajectory spaces : Accelerated thermodynamics and kinetic calculations
////本方法中的A,B两类原子将被称作蛋白,水。
//
//struct SITS_INFORMATION
//{
//	int atom_numbers;//体系的所有原子数目,跟随MD_INFORMATION中的值
//	int protein_atom_numbers;//体系的蛋白原子数目,人为设定
//	int water_atom_numbers;//体系的水原子数目
//
//	NEIGHBOR_LIST *protein_water_neighbor;//蛋白原子周围的水邻居表，长度为protein_atom_numbers
//	int *protein_water_numbers;//蛋白周围水邻居的原子个数，长度为protein_atom_numbers
//	//蛋白-蛋白，水-水近邻表直接用普通MD中声明的近邻表
//
//	VECTOR *protein_water_frc;//蛋白-水原子间的frc，长度为atom_numbers，前protein_atom_numbers个为蛋白原子受到水的frc分量，后部分为水受到蛋白原子的frc分量
//	//蛋白-蛋白，水-水作用的frc直接用普通MD中生命的frc
//
//	float *protein_water_energy;//蛋白与水原子间的作用能量，长度为protein_atom_numbers
//	//假设蛋白-水之间不存在键的作用而只有LJ势能和CF势能，因此蛋白-蛋白能量，水-水能量可按照未区分的情况直接存于一个列表中
//	//由于不区分某个部分内部的具体的能量分量，因此可以直接把bond angle和dihedral能量加到参与的任意一个原子身上
//	float *atom_energy;//长度为atom_numbers，包含蛋白-蛋白能量，水-水能量
//
//	//长度均为1的用于存储总的分块能量的指针
//	float *sum_of_water_water_energy;
//	float *sum_of_protein_protein_energy;
//	float *sum_of_protein_water_energy;
//
//
//	int k_numbers;
//	float *beta_k;
//	float *NkExpBetakU;
//	float *Nk;
//	float *sum_a;
//	float *sum_b;
//	float *factor;
//
//	//迭代使用，XYJ
//	//控制变量
//	
//	int record_count;
//	int record_interval;
//	int update_interval;
//	int reset = 1;
//
//	//计算变量
//	//xyj的cpp变量名-ylj的F90变量名-文献对应
//	//ene_recorded - vshift - ene
//	//gf - gf - log( n_k * exp(-beta_k * ene) )
//	//gfsum - gfsum - log( Sum_(k=1)^N ( log( n_k * exp(-beta_k * ene) ) ) )
//	//log_weight - rb - log of the weighting function
//	//log_mk_inverse - ratio - log(m_k^-1)
//	//log_norm_old - normlold - W(j-1)
//	//log_norm - norml - W(j)
//	//log_pk - rbfb - log(p_k)
//	//log_nk_inverse - pratio - log(n_k^-1)
//	//log_nk - fb - log(n_k)
//
//	float *ene_recorded;
//	float *gf;
//	float *gfsum;
//	float *log_weight;
//	float *log_mk_inverse;
//	float *log_norm_old;
//	float *log_norm;
//	float *log_pk;
//	float *log_nk_inverse;
//	float *log_nk;
//
//
//};
//
////记录时更新函数
//
//__global__ void SITS_Record_Ene(float *ene_record, const float *pw_ene, const float *pp_ene, const float pe_a, const float pe_b);
//
//__global__ void SITS_Update_gf(const int kn, float *gf,
//	const float *ene_record, const float *log_nk, const float *beta_k);
//
//__global__ void SITS_Update_gfsum(const int kn, float *gfsum, const float *gf);
//
//__global__ void SITS_Update_log_pk(const int kn, float *log_pk,
//	const float *gf, const float *gfsum, const int reset);
//
////迭代时更新函数
//
//__global__ void SITS_Update_log_mk_inverse(const int kn,
//	float *log_weight, float *log_mk_inverse, float *log_norm_old,
//	float *log_norm, const float *log_pk, const float *log_nk);
//
//__global__ void SITS_Update_log_nk_inverse(const int kn,
//	float *log_nk_inverse, const float *log_mk_inverse);
//
//__global__ void SITS_Update_nk(const int kn,
//	float *log_nk, float *nk, const float *log_nk_inverse);
//
//
//
//#endif //SITS_COMMON_CUH(SITS_common.cuh)