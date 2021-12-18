//#ifndef SITS_MD_FORCE_CUH
//#define SITS_MD_FORCE_CUH
//#include "SITS_common.cuh"
////����ĺ������㷽������ͨMD�ж�Ӧ�ĸ�����û��̫���𣬣�ע��Direct����ʹ����PME�����е�ʽ�ӣ���SITS�д��ɣ�
////ÿ������ÿ�μ��㶼Ҫ�õ������������������Ǽӵ�ԭ�����ϣ����ڲ�����������������˼ӵ�ÿ����������ĵ�һ��ԭ���ϣ�
////�������һ����Ҫ�ý��ڱ�����LJ��CF�⣬����������ֱ�Ӷ�frc�б��ԭ�������б����ã����Զ�����Ϊ����-���ײ��ֺ�ˮ-ˮ���֣�
////���һ������Ҫ�������Σ�һ���Ǽ��㵰��-���ף�ˮ-ˮ�໥���ú��������������������죻�ڶ��μ���ֻ���㵽protein_atom_numbers��
////��������������ֱ����һ���µ�frc�У���¼ˮ�Ե��׵����Լ����׶�ˮ��������energy list�У�����Ϊprotein_atom_numbers����¼�˵���-ˮ���໥���ã�
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
////xyj�޸Ĺ��������pe_a��pe_b�������Ա仯������fb_bias
//__global__ void SITS_Enhanced_Force
//(const int atom_numbers, const int protein_atom_numbers,
//VECTOR *md_frc, const VECTOR *pw_frc,
//const float *pp_ene, const float *pw_ene,
//const int k_numbers, float *nkexpbetaku,
//const float *beta_k, const float *n_k,
//float *sum_a, float *sum_b, float *factor,
//const float beta_0, const float pe_a, const float pe_b, const float fb_bias);
//#endif //SITS_MD_FORCE_CUH(SITS_MD_force.cuh)