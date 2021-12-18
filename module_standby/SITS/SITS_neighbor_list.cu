//#include "SITS_neighbor_list.cuh"
//__global__ void SITS_Put_Atom_In_Grid_Bucket(const int atom_start,const int atom_end, const int *atom_in_grid_serial, GRID_BUCKET *bucket, int *atom_numbers_in_grid_bucket)
//{
//	int atom_i = blockDim.x*blockIdx.x + threadIdx.x + atom_start;
//	if (atom_i < atom_end)
//	{
//		int grid_serial = atom_in_grid_serial[atom_i];
//		int a = atom_numbers_in_grid_bucket[grid_serial];
//		atomicCAS(&bucket[grid_serial].atom_serial[a], -1, atom_i);
//		if (bucket[grid_serial].atom_serial[a] != atom_i)
//		{
//			while (true)
//			{
//				a = a + 1;
//				atomicCAS(&bucket[grid_serial].atom_serial[a], -1, atom_i);
//				if (bucket[grid_serial].atom_serial[a] == atom_i)
//				{
//					atomicAdd(&atom_numbers_in_grid_bucket[grid_serial], 1);
//					break;
//				}
//			}
//		}
//		else
//		{
//			atomicAdd(&atom_numbers_in_grid_bucket[grid_serial], 1);
//		}
//	}
//}
//__global__ void SITS_Find_atom_neighbors(
//	const int atom_start,const int atom_end, const UNSIGNED_INT_VECTOR *uint_crd, const VECTOR uint_dr_to_dr_cof,
//	const int *atom_in_grid_serial, const GRID_POINTER *gpointer, const GRID_BUCKET *bucket, const int *atom_numbers_in_grid_bucket,
//	NEIGHBOR_LIST *nl, int *atomnumbers_in_nl, const float cutoff_skin_square)
//{
//	int atom_i = blockDim.x*blockIdx.x + threadIdx.x + atom_start;
//	if (atom_i < atom_end)
//	{
//		int grid_serial = atom_in_grid_serial[atom_i];
//		int grid_serial2;
//		int atom_numbers_in_nl_lin = 0;
//		int atom_j;
//		int int_x;
//		int int_y;
//		int int_z;
//		VECTOR dr;
//		float dr2;
//		for (int grid_cycle = 0; grid_cycle < 125; grid_cycle = grid_cycle + 1)
//		{
//			grid_serial2 = gpointer[grid_serial].grid_serial[grid_cycle];
//			for (int i = 0; i < atom_numbers_in_grid_bucket[grid_serial2]; i = i + 1)
//			{
//				atom_j = bucket[grid_serial2].atom_serial[i];
//				if (atom_j > atom_i)
//				{
//					int_x = uint_crd[atom_j].uint_x - uint_crd[atom_i].uint_x;
//					int_y = uint_crd[atom_j].uint_y - uint_crd[atom_i].uint_y;
//					int_z = uint_crd[atom_j].uint_z - uint_crd[atom_i].uint_z;
//					dr.x = uint_dr_to_dr_cof.x*int_x;
//					dr.y = uint_dr_to_dr_cof.y*int_y;
//					dr.z = uint_dr_to_dr_cof.z*int_z;
//					dr2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
//
//					if (dr2 < cutoff_skin_square)
//					{
//						nl[atom_i].atom_serial[atom_numbers_in_nl_lin] = atom_j;
//						atom_numbers_in_nl_lin = atom_numbers_in_nl_lin + 1;
//					}
//				}
//			}
//		}//124 grid cycle
//		atomnumbers_in_nl[atom_i] = atom_numbers_in_nl_lin;
//	}
//}
//__global__ void SITS_Find_protein_water_neighbors(
//	const int atom_end, const UNSIGNED_INT_VECTOR *uint_crd, const VECTOR uint_dr_to_dr_cof,
//	const int *atom_in_grid_serial, const GRID_POINTER *gpointer, const GRID_BUCKET *bucket, const int *atom_numbers_in_grid_bucket,
//	NEIGHBOR_LIST *nl, int *atomnumbers_in_nl, const float cutoff_skin_square)
//{
//	int atom_i = blockDim.x*blockIdx.x + threadIdx.x;
//	if (atom_i < atom_end)
//	{
//		int grid_serial = atom_in_grid_serial[atom_i];
//		int grid_serial2;
//		int atom_numbers_in_nl_lin = 0;
//		int atom_j;
//		int int_x;
//		int int_y;
//		int int_z;
//		VECTOR dr;
//		float dr2;
//		for (int grid_cycle = 0; grid_cycle < 125; grid_cycle = grid_cycle + 1)
//		{
//			grid_serial2 = gpointer[grid_serial].grid_serial[grid_cycle];
//			for (int i = 0; i < atom_numbers_in_grid_bucket[grid_serial2]; i = i + 1)
//			{
//				atom_j = bucket[grid_serial2].atom_serial[i];
//				if (atom_j >=atom_end)
//				{
//					int_x = uint_crd[atom_j].uint_x - uint_crd[atom_i].uint_x;
//					int_y = uint_crd[atom_j].uint_y - uint_crd[atom_i].uint_y;
//					int_z = uint_crd[atom_j].uint_z - uint_crd[atom_i].uint_z;
//					dr.x = uint_dr_to_dr_cof.x*int_x;
//					dr.y = uint_dr_to_dr_cof.y*int_y;
//					dr.z = uint_dr_to_dr_cof.z*int_z;
//					dr2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
//
//					if (dr2 < cutoff_skin_square)
//					{
//						nl[atom_i].atom_serial[atom_numbers_in_nl_lin] = atom_j;
//						atom_numbers_in_nl_lin = atom_numbers_in_nl_lin + 1;
//					}
//				}
//			}
//		}//124 grid cycle
//		atomnumbers_in_nl[atom_i] = atom_numbers_in_nl_lin;
//	}
//}
//
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
//int *excluded_list_start, int * excluded_list, int * excluded_numbers)
//{
//	if (refresh_sign[0] == 1)
//	{
//		int p_atom_numbers=atom_end;
//		int w_atom_numbers=atom_numbers-atom_end;
//		/**********************************************************************************************************/
//		//与普通MD一样的更新过程，先进行向量平移，做周期性映射，根据应放入的桶子对原子编号，重新平移坐标，更新老坐标，得到整数坐标
//		VECTOR trans_vec = { -half_skin_square, -half_skin_square, -half_skin_square };
//		Vector_Translation << <ceilf((float)atom_numbers / thread), thread >> >(atom_numbers, crd, trans_vec);
//		Crd_Periodic_Map << <ceilf((float)atom_numbers / thread), thread >> >(atom_numbers, crd, box_length);
//		Find_Atom_In_Grid_Serial << <ceilf((float)atom_numbers / thread), thread >> >
//			(atom_numbers, grid_info.grid_length_inverse, crd, grid_info.grid_N, grid_info.Nxy, atom_in_grid_serial);
//		trans_vec.x = -trans_vec.x;
//		trans_vec.y = -trans_vec.y;
//		trans_vec.z = -trans_vec.z;
//		Vector_Translation << <ceilf((float)atom_numbers / thread), thread >> >(atom_numbers, crd, trans_vec);
//		Copy_List << <ceilf((float)3.*atom_numbers / thread), thread >> >
//			(3 * atom_numbers, (float*)crd, (float*)old_crd);
//		Crd_To_Uint_Crd << <ceilf((float)atom_numbers / thread), thread >> >
//			(atom_numbers, crd_to_uint_crd_cof, crd, uint_crd);
//		/**********************************************************************************************************/
//
//		//清空桶子中的原子
//		Clear_Grid_Bucket << <ceilf((float)grid_info.grid_numbers / thread), thread >> >
//			(grid_info.grid_numbers, atom_numbers_in_grid_bucket, bucket);
//		//向桶子中放入蛋白原子
//		SITS_Put_Atom_In_Grid_Bucket << <ceilf((float)atom_end / thread), thread >> >
//			(0,atom_end, atom_in_grid_serial, bucket, atom_numbers_in_grid_bucket);
//		//构建蛋白-蛋白近邻表
//		SITS_Find_atom_neighbors << <ceilf((float)atom_end / thread), thread >> >
//			(0,atom_end, uint_crd, uint_dr_to_dr_cof,
//			atom_in_grid_serial, gpointer, bucket, atom_numbers_in_grid_bucket,
//			d_nl, d_atomnumber_in_nl, 100);//144
//
//		//清空桶子中的原子
//		/*Clear_Grid_Bucket << <ceilf((float)grid_info.grid_numbers / thread), thread >> >
//			(grid_info.grid_numbers, atom_numbers_in_grid_bucket, bucket);*/
//		//向桶子中放入水原子
//		SITS_Put_Atom_In_Grid_Bucket << <ceilf((float)w_atom_numbers / thread), thread >> >
//			(atom_end, atom_numbers, atom_in_grid_serial, bucket, atom_numbers_in_grid_bucket);
//		//构建水-水近邻表
//		SITS_Find_atom_neighbors << <ceilf((float)w_atom_numbers / thread), thread >> >
//			(atom_end, atom_numbers, uint_crd, uint_dr_to_dr_cof,
//			atom_in_grid_serial, gpointer, bucket, atom_numbers_in_grid_bucket,
//			d_nl, d_atomnumber_in_nl, 100);//144
//
//		//进行蛋白-蛋白，水-水近邻表中的原子剔除
//		Delete_Excluded_Atoms_Serial_In_Neighbor_List << <ceilf((float)atom_numbers / thread), thread >> >
//			(atom_numbers, d_nl, d_atomnumber_in_nl, excluded_list_start, excluded_list, excluded_numbers);
//
//
//		//向桶子中再放入蛋白原子
//		/*SITS_Put_Atom_In_Grid_Bucket << <ceilf((float)atom_end / thread), thread >> >
//			(0, atom_end, atom_in_grid_serial, bucket, atom_numbers_in_grid_bucket);*/
//		//构建蛋白-水近邻表。注意，此时构建的近邻表是非全原子的近邻表而是对蛋白原子一个专用的近邻表（用于存储蛋白-水交叉项）
//		SITS_Find_protein_water_neighbors << <ceilf((float)atom_end / thread), thread >> >
//			(atom_end, uint_crd, uint_dr_to_dr_cof,
//			atom_in_grid_serial, gpointer, bucket, atom_numbers_in_grid_bucket,
//			d_nl_p, d_atomnumber_in_nl_p, 100);//144
//		//假设蛋白和水原子间没有剔除关系，因此不对这个近邻表进行剔除操作。
//
//		refresh_sign[0] = 0;
//	}
//}