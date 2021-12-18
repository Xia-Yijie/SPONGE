//#ifndef SITS_MD_FORCE_CUH
//#define SITS_MD_FORCE_CUH
//#include "SITS_common.cuh"
////下面的函数计算方法与普通MD中对应的个函数没有太大差别，（注意Direct部分使用了PME方法中的式子，在SITS中存疑）
////每个函数每次计算都要得到力和能量，能量总是加到原子身上（由于不区分能量分量，因此加到每个能量计算的第一个原子上）
////除了最后一个需要用近邻表计算的LJ和CF外，其他函数均直接对frc列表和原子能量列表作用（其自动区分为蛋白-蛋白部分和水-水部分）
////最后一个函数要调用两次，一次是计算蛋白-蛋白，水-水相互作用和能量，与其他函数无异；第二次计算只计算到protein_atom_numbers，
////算出的力和能量分别存于一个新的frc中（记录水对蛋白的力以及蛋白对水的力）和energy list中（长度为protein_atom_numbers，记录了蛋白-水的相互作用）
//
//__global__ void SITS_Angle_Force_Energy
//(const int angle_numbers,
//const UNSIGNED_INT_VECTOR *uint_crd, const VECTOR scaler,
//const int *atom_a, const int *atom_b, const int *atom_c,
//const float *angle_k, const float *angle_theta0,
//VECTOR *frc, float *energy);
//
//__global__ void SITS_Bond_Force_Energy
//(const int bond_numbers, const UNSIGNED_INT_VECTOR *uint_crd, const VECTOR scaler,
//const int *atom_a, const int *atom_b, const float *bond_k, const float *bond_r0, VECTOR *frc, float *energy);
//
//__global__ void SITS_Dihedral_Force_Energy
//(const int dihedral_numbers, const UNSIGNED_INT_VECTOR *uint_crd, const VECTOR scaler,
//const int *atom_a, const int *atom_b, const int *atom_c, const int *atom_d, const int *ipn, const float *pk, const float *gamc,
//const float *gams, const float *pn, VECTOR *frc, float *energy);
//
//__global__ void SITS_Dihedral_14_LJ_Force_With_Direct_CF_Energy
//(const int dihedral_14_numbers, const UINT_VECTOR_LJ_TYPE *uint_crd, const VECTOR boxlength,
//const int *a_14, const int *b_14, const float *lj_scale_factor, const float *cf_scale_factor, const float *LJ_type_A, const float *LJ_type_B, VECTOR *frc, float *energy);
//
//__global__ void SITS_LJ_Force_With_Direct_CF_Energy
//(const int atom_numbers, const NEIGHBOR_LIST *nl, const int *atomnumbers_in_nl,
//const UINT_VECTOR_LJ_TYPE *uint_crd, const VECTOR boxlength,
//const float *LJ_type_A, const float *LJ_type_B, const float cutoff,
//VECTOR *frc, const float pme_beta, const float sqrt_pi, float *energy);
//
//
////xyj修改过，添加了pe_a、pe_b两个线性变化变量和fb_bias
//__global__ void SITS_Enhanced_Force
//(const int atom_numbers, const int protein_atom_numbers,
//VECTOR *md_frc, const VECTOR *pw_frc,
//const float *pp_ene, const float *pw_ene,
//const int k_numbers, float *nkexpbetaku,
//const float *beta_k, const float *n_k,
//float *sum_a, float *sum_b, float *factor,
//const float beta_0, const float pe_a, const float pe_b, const float fb_bias);
//#endif //SITS_MD_FORCE_CUH(SITS_MD_force.cuh)