//#ifndef SITS_NEIGHBOR_LIST
//#define SITS_NEIGHBOR_LIST
//#include "SITS_common.cuh"
////功能与原始的Put_Atom_In_Grid_Bucket函数类似
////差别在于只放入原子编号为[atom_start,atom_end)间的原子进入筒子
////目的在于减少构建蛋白-蛋白、水-水近邻表过程中的计算量
////在清空的Bucket基础上先放入蛋白，构建蛋白-蛋白近邻表；再清空桶子，放入水，构建水-水近邻表；不清空桶子，再放入蛋白，用于构建蛋白-水近邻表
//__global__ void SITS_Put_Atom_In_Grid_Bucket
//(const int atom_start, const int atom_end, const int *atom_in_grid_serial, GRID_BUCKET *bucket, int *atom_numbers_in_grid_bucket);
//
////用于构建蛋白-蛋白、水-水近邻表，原子编号为[atom_start,atom_end)间的原子
//__global__ void SITS_Find_atom_neighbors
//(const int atom_start, const int atom_end, const UNSIGNED_INT_VECTOR *uint_crd, const VECTOR uint_dr_to_dr_cof,
//const int *atom_in_grid_serial, const GRID_POINTER *gpointer, const GRID_BUCKET *bucket, const int *atom_numbers_in_grid_bucket,
//NEIGHBOR_LIST *nl, int *atomnumbers_in_nl, const float cutoff_skin_square);
//
////用于构建蛋白-水近邻表，默认蛋白的原子编号为[0,atom_end)区间，而水原子的编号为[atom_end,atom_numbers)区间
////因此，该近邻表的长度应只有atom_end个，不要混用全原子的近邻表
//__global__ void SITS_Find_protein_water_neighbors
//(const int atom_end, const UNSIGNED_INT_VECTOR *uint_crd, const VECTOR uint_dr_to_dr_cof,
//const int *atom_in_grid_serial, const GRID_POINTER *gpointer, const GRID_BUCKET *bucket, const int *atom_numbers_in_grid_bucket,
//NEIGHBOR_LIST *nl, int *atomnumbers_in_nl, const float cutoff_skin_square);
//
////该函数的具体解释，请看函数的定义
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