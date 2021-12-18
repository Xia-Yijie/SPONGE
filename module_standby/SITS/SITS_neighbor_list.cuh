//#ifndef SITS_NEIGHBOR_LIST
//#define SITS_NEIGHBOR_LIST
//#include "SITS_common.cuh"
////������ԭʼ��Put_Atom_In_Grid_Bucket��������
////�������ֻ����ԭ�ӱ��Ϊ[atom_start,atom_end)���ԭ�ӽ���Ͳ��
////Ŀ�����ڼ��ٹ�������-���ס�ˮ-ˮ���ڱ�����еļ�����
////����յ�Bucket�������ȷ��뵰�ף���������-���׽��ڱ������Ͱ�ӣ�����ˮ������ˮ-ˮ���ڱ������Ͱ�ӣ��ٷ��뵰�ף����ڹ�������-ˮ���ڱ�
//__global__ void SITS_Put_Atom_In_Grid_Bucket
//(const int atom_start, const int atom_end, const int *atom_in_grid_serial, GRID_BUCKET *bucket, int *atom_numbers_in_grid_bucket);
//
////���ڹ�������-���ס�ˮ-ˮ���ڱ�ԭ�ӱ��Ϊ[atom_start,atom_end)���ԭ��
//__global__ void SITS_Find_atom_neighbors
//(const int atom_start, const int atom_end, const UNSIGNED_INT_VECTOR *uint_crd, const VECTOR uint_dr_to_dr_cof,
//const int *atom_in_grid_serial, const GRID_POINTER *gpointer, const GRID_BUCKET *bucket, const int *atom_numbers_in_grid_bucket,
//NEIGHBOR_LIST *nl, int *atomnumbers_in_nl, const float cutoff_skin_square);
//
////���ڹ�������-ˮ���ڱ�Ĭ�ϵ��׵�ԭ�ӱ��Ϊ[0,atom_end)���䣬��ˮԭ�ӵı��Ϊ[atom_end,atom_numbers)����
////��ˣ��ý��ڱ�ĳ���Ӧֻ��atom_end������Ҫ����ȫԭ�ӵĽ��ڱ�
//__global__ void SITS_Find_protein_water_neighbors
//(const int atom_end, const UNSIGNED_INT_VECTOR *uint_crd, const VECTOR uint_dr_to_dr_cof,
//const int *atom_in_grid_serial, const GRID_POINTER *gpointer, const GRID_BUCKET *bucket, const int *atom_numbers_in_grid_bucket,
//NEIGHBOR_LIST *nl, int *atomnumbers_in_nl, const float cutoff_skin_square);
//
////�ú����ľ�����ͣ��뿴�����Ķ���
//__global__ void SITS_Refresh_Neighbor_List
//(int *refresh_sign, const int thread,
//const int atom_numbers, const int atom_end,
//VECTOR *crd, VECTOR *old_crd, UNSIGNED_INT_VECTOR *uint_crd,
//const VECTOR crd_to_uint_crd_cof, const VECTOR uint_dr_to_dr_cof,
//int *atom_in_grid_serial,
//const float half_skin_square, const VECTOR box_length,
//const GRID_INFORMATION grid_info, const GRID_POINTER *gpointer,
//GRID_BUCKET *bucket, int *atom_numbers_in_grid_bucket,
//
//NEIGHBOR_LIST *d_nl, int *d_atomnumber_in_nl,
//NEIGHBOR_LIST *d_nl_p, int *d_atomnumber_in_nl_p,
//
//int *excluded_list_start, int * excluded_list, int * excluded_numbers);
//
//#endif //SITS_NEIGHBOR_LIST(SITS_neighbor_list.cuh)