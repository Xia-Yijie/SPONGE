//#include "SITS_common.cuh"
//#include "SITS_MD_force.cuh"
//#include "SITS_neighbor_list.cuh"
//
//__global__ void Add_List(const int element_numbers, float *list1, const float *list2)
//{
//	for (int i = threadIdx.x; i < element_numbers; i = i + blockDim.x)
//	{
//		list1[i] = list1[i] + list2[i];
//	}
//}
//__global__ void Sum_Of_List(const int start,const int end, const float* list, float *sum)
//{
//	if (threadIdx.x == 0)
//	{
//		sum[0] = 0.;
//	}
//	__syncthreads();
//	float lin = 0.;
//	for (int i = threadIdx.x + start; i < end; i = i + blockDim.x)
//	{
//		lin = lin + list[i];
//	}
//	atomicAdd(sum, lin);
//}
//#include <Windows.h>//不可少 
//#include <mmsystem.h>//不可少
//int main_SITS(AMBER_SITS_INFORMATION amber_sits_info)
//{
//	AMBER_INFORMATION amber_info = amber_sits_info.amber_info;
//	MD_INFORMATION md_info;//声明一个MD基本信息
//	LENNARD_JONES_INFOMATION LJ_info;
//	NON_BOND_INFORMATION nb_info;
//	GRID_INFORMATION grid_info;
//	CONSTANT const_parameter;
//	md_info.dt = 0.020455 * amber_info.dt;
//	BOND_INFORMATION bond_info;
//	Import_Bond_Information_From_AMBERFILE(amber_info.Parm_File_Name,
//		&bond_info.bond_numbers, &bond_info.atom_a, &bond_info.atom_b, &bond_info.k, &bond_info.r0);
//	ANGLE_INFORMATION angle_info;
//	Import_Angle_Information_From_AMBERFILE(amber_info.Parm_File_Name,
//		&angle_info.angle_numbers, &angle_info.atom_a, &angle_info.atom_b, &angle_info.atom_c,
//		&angle_info.k, &angle_info.theta0);
//	DIHEDRAL_INFORMATION dihedral_info;
//	Import_Dihedral_Information_From_AMBERFILE(amber_info.Parm_File_Name,
//		&dihedral_info.dihedral_numbers, &dihedral_info.dihedral_14_numbers,
//		&dihedral_info.atom_a, &dihedral_info.atom_b, &dihedral_info.atom_c, &dihedral_info.atom_d,
//		&dihedral_info.pk, &dihedral_info.gamc, &dihedral_info.gams, &dihedral_info.pn, &dihedral_info.ipn,
//		&dihedral_info.atom_a_14, &dihedral_info.atom_b_14, &dihedral_info.lj_scale_factor, &dihedral_info.cf_scale_factor);
//
//
//	//从parm7中读入基本信息
//	Import_Basic_System_Information_From_AMBERFILE(amber_info.Parm_File_Name,
//		&md_info.atom_numbers,
//		&md_info.residue_numbers, &md_info.res_start, &md_info.res_end,
//		&md_info.d_mass, &md_info.d_charge);
//
//	Malloc_Safely((void**)&md_info.h_mass, sizeof(float)*md_info.atom_numbers);
//	cudaMemcpy(md_info.h_mass, md_info.d_mass, sizeof(float)*md_info.atom_numbers, cudaMemcpyDeviceToHost);
//	LJ_info.Initial_From_AMBER_Parm(amber_info.Parm_File_Name);
//	Import_Excluded_List_From_AMBERFILE(amber_info.Parm_File_Name, &nb_info.excluded_list_start, &nb_info.excluded_list, &nb_info.excluded_numbers);
//	Malloc_Safely((void**)&md_info.coordinate, sizeof(VECTOR)*md_info.atom_numbers);
//	Malloc_Safely((void**)&md_info.velocity, sizeof(VECTOR)*md_info.atom_numbers);
//	Malloc_Safely((void**)&md_info.force, sizeof(VECTOR)*md_info.atom_numbers);
//	//从rst7中读入基本信息
//	Import_Information_From_Rst7(amber_info.Crd_File_Name,
//		&md_info.atom_numbers, &md_info.simulation_start_time, &md_info.crd, &md_info.vel, &md_info.box_length, amber_info.irest);
//	//初始化另外一些基本信息
//	Cuda_Malloc_Safely((void**)&md_info.old_crd, sizeof(VECTOR)*md_info.atom_numbers);
//	Cuda_Malloc_Safely((void**)&md_info.acc, sizeof(VECTOR)*md_info.atom_numbers);
//	Cuda_Malloc_Safely((void**)&md_info.frc, sizeof(VECTOR)*md_info.atom_numbers);
//	Cuda_Malloc_Safely((void**)&md_info.d_mass_inverse, sizeof(float)*md_info.atom_numbers);
//	Reset_List << <ceilf((float)3.*md_info.atom_numbers / 32), 32 >> >
//		(3 * md_info.atom_numbers, (float*)md_info.acc, 0.);
//	Reset_List << <ceilf((float)3.*md_info.atom_numbers / 32), 32 >> >
//		(3 * md_info.atom_numbers, (float*)md_info.frc, 0.);
//	Inverse_List_Element << <ceilf((float)md_info.atom_numbers / 32), 32 >> >
//		(md_info.atom_numbers, md_info.d_mass, md_info.d_mass_inverse);
//
//	//初始化近邻表和grid等信息
//	Initial_Neighbor_List(md_info.atom_numbers, &nb_info.h_nl, &nb_info.d_nl, &nb_info.h_atomnumbers_in_nl, &nb_info.d_atomnumbers_in_nl, amber_info.max_neibor_numbers);
//	printf("Initial Neighbor Successed!\n");
//
//	nb_info.atom_numbers = &md_info.atom_numbers; //设定一些模拟信息
//	nb_info.box_length = &md_info.box_length;
//	nb_info.half_box_length = &md_info.half_box_length;
//	nb_info.cutoff = amber_info.cut;
//	nb_info.skin = amber_info.skin;
//	Cuda_Malloc_Safely((void**)&nb_info.is_need_refresh_neighbor_list, sizeof(int));
//	Reset_List << <1, 1 >> >(1, nb_info.is_need_refresh_neighbor_list, 0);
//
//	Cuda_Malloc_Safely((void**)&grid_info.atom_in_grid_serial, sizeof(int)*md_info.atom_numbers);
//	Initial_Neighbor_Grid(
//		&grid_info.gpointer, &grid_info.bucket, &grid_info.atom_numbers_in_grid_bucket,
//		&nb_info, &grid_info,
//		amber_info.max_atom_in_grid_numbers);
//	printf("Initial Neighbor_Grid Successed!\n");
//
//	md_info.crd_to_uint_crd_cof = { const_parameter.uint_max_float / md_info.box_length.x, const_parameter.uint_max_float / md_info.box_length.y, const_parameter.uint_max_float / md_info.box_length.z };
//	md_info.uint_dr_to_dr_cof = { 1. / md_info.crd_to_uint_crd_cof.x, 1. / md_info.crd_to_uint_crd_cof.y, 1. / md_info.crd_to_uint_crd_cof.z };
//	md_info.uint_crd_to_grid_serial = { md_info.uint_dr_to_dr_cof.x*grid_info.grid_length_inverse.x,
//		md_info.uint_dr_to_dr_cof.y*grid_info.grid_length_inverse.y,
//		md_info.uint_dr_to_dr_cof.z*grid_info.grid_length_inverse.z };
//	Cuda_Malloc_Safely((void **)&md_info.uint_crd, sizeof(UNSIGNED_INT_VECTOR)*md_info.atom_numbers);
//	Cuda_Malloc_Safely((void**)&md_info.uint_crd_with_LJ, sizeof(UINT_VECTOR_LJ_TYPE)*md_info.atom_numbers);
//	printf("Initial Uint_Crd Successed!\n");
//
//	VECTOR trans_vec = { amber_info.skin, amber_info.skin, amber_info.skin };
//	VECTOR trans_vec_minus = -1.0 * trans_vec;
//	Clear_Grid_Bucket << <ceilf((float)grid_info.grid_numbers / 32), 32 >> >
//		(grid_info.grid_numbers, grid_info.atom_numbers_in_grid_bucket, grid_info.bucket);
//	Crd_Periodic_Map << <ceilf((float)md_info.atom_numbers / 32), 32 >> >(md_info.atom_numbers, md_info.crd, md_info.box_length);
//
//	Find_Atom_In_Grid_Serial << <ceilf((float)md_info.atom_numbers / 32), 32 >> >
//		(md_info.atom_numbers, grid_info.grid_length_inverse, md_info.crd, grid_info.grid_N, grid_info.Nxy, grid_info.atom_in_grid_serial);
//
//	Vector_Translation << <ceilf((float)md_info.atom_numbers / 32), 32 >> >(md_info.atom_numbers, md_info.crd, trans_vec);
//	Copy_List << <ceilf((float)3.*md_info.atom_numbers / 32), 32 >> >
//		(3 * md_info.atom_numbers, (float*)md_info.crd, (float*)md_info.old_crd);
//
//	Put_Atom_In_Grid_Bucket << <ceilf((float)md_info.atom_numbers / 32), 32 >> >
//		(md_info.atom_numbers, grid_info.atom_in_grid_serial, grid_info.bucket, grid_info.atom_numbers_in_grid_bucket);
//
//	Crd_To_Uint_Crd << <ceilf((float)md_info.atom_numbers / 32), 32 >> >
//		(md_info.atom_numbers, 0.5*md_info.crd_to_uint_crd_cof, md_info.crd, md_info.uint_crd);
//	Find_atom_neighbors << <ceilf((float)md_info.atom_numbers / 32), 32 >> >
//		(md_info.atom_numbers, md_info.uint_crd, md_info.uint_dr_to_dr_cof,
//		grid_info.atom_in_grid_serial, grid_info.gpointer, grid_info.bucket, grid_info.atom_numbers_in_grid_bucket,
//		nb_info.d_nl, nb_info.d_atomnumbers_in_nl, nb_info.cutoff_with_skin_square);
//	Delete_Excluded_Atoms_Serial_In_Neighbor_List << <ceilf((float)md_info.atom_numbers / 32), 32 >> >
//		(md_info.atom_numbers, nb_info.d_nl, nb_info.d_atomnumbers_in_nl, nb_info.excluded_list_start, nb_info.excluded_list, nb_info.excluded_numbers);
//
//	//朗之万热浴
//	//Malloc_Safely((void**)&md_info.h_mass, sizeof(float)*md_info.atom_numbers);
//	//cudaMemcpy(md_info.h_mass, md_info.d_mass, sizeof(float)*md_info.atom_numbers, cudaMemcpyDeviceToHost);
//	LIUJIAN_MD_INFORMATION liujian_info;
//	liujian_info.atom_numbers = md_info.atom_numbers;
//	liujian_info.dt = md_info.dt;
//	liujian_info.gamma_ln = amber_info.gamma_ln;
//	liujian_info.target_temperature = amber_info.temp0;
//	Initial_LIUJIAN_MD(&liujian_info, md_info.h_mass);
//	Setup_Rand_Normal_Kernel << <ceilf((float)liujian_info.float4_numbers / 32.), 32 >> >
//		(liujian_info.float4_numbers, liujian_info.rand_state, 0);
//
//
//	//初始化PME
//	Particle_Mesh_Ewald pme_info;
//	pme_info.Control_From_AMBER_Mdin(amber_info.In_File_Name);
//	pme_info.Initial_PME_System(md_info.atom_numbers,md_info.box_length);
//	
//
//	//能量相关
//	Cuda_Malloc_Safely((void**)&md_info.res_ek_energy, sizeof(float)*md_info.residue_numbers);
//	Cuda_Malloc_Safely((void**)&md_info.sigma_of_res_ek, sizeof(float)* 1);
//	Cuda_Malloc_Safely((void**)&bond_info.bond_ene, sizeof(float)*bond_info.bond_numbers);
//	Cuda_Malloc_Safely((void**)&bond_info.sigme_of_bond_ene, sizeof(float)* 1);
//	Cuda_Malloc_Safely((void**)&angle_info.energy, sizeof(float)*angle_info.angle_numbers);
//	Cuda_Malloc_Safely((void**)&angle_info.sigma_energy, sizeof(float)* 1);
//	
//	Cuda_Malloc_Safely((void**)&dihedral_info.dihedral_energy, sizeof(float)*dihedral_info.dihedral_numbers);
//	Cuda_Malloc_Safely((void**)&dihedral_info.dihedral_14_energy, sizeof(float)*dihedral_info.dihedral_14_numbers);
//	Cuda_Malloc_Safely((void**)&dihedral_info.sigma_energy, sizeof(float)* 1);
//	Cuda_Malloc_Safely((void**)&dihedral_info.sigma_14_LJ_energy, sizeof(float)* 1);
//	Cuda_Malloc_Safely((void**)&dihedral_info.sigma_14_CF_energy, sizeof(float)* 1);
//
//	Copy_Crd_To_New_Crd_Start << <ceilf((float)md_info.atom_numbers / 32), 32 >> >
//		(md_info.atom_numbers, md_info.uint_crd, md_info.uint_crd_with_LJ, LJ_info.d_atom_LJ_type, md_info.d_charge);
//
//
//
//
//	//SITS
//	/****************************************************************************************************/
//	SITS_INFORMATION sits_info;
//	sits_info.atom_numbers = md_info.atom_numbers;
//	sits_info.protein_atom_numbers = amber_sits_info.selected_atom_numbers;
//	sits_info.water_atom_numbers = sits_info.atom_numbers - sits_info.protein_atom_numbers;
//
//	NEIGHBOR_LIST *h_lin;
//	Malloc_Safely((void**)&h_lin, sizeof(NEIGHBOR_LIST)*sits_info.protein_atom_numbers);
//	for (int i = 0; i < sits_info.protein_atom_numbers; i = i + 1)
//	{
//		Cuda_Malloc_Safely((void**)&h_lin[i].atom_serial, sizeof(int)* amber_info.max_neibor_numbers);
//	}
//	Cuda_Malloc_Safely((void**)&sits_info.protein_water_neighbor, sizeof(NEIGHBOR_LIST)*sits_info.protein_atom_numbers);
//	cudaMemcpy(sits_info.protein_water_neighbor, h_lin, sizeof(NEIGHBOR_LIST)*sits_info.protein_atom_numbers,cudaMemcpyHostToDevice);
//	free(h_lin);
//	Cuda_Malloc_Safely((void**)&sits_info.protein_water_numbers, sizeof(int)*sits_info.protein_atom_numbers);
//
//	Cuda_Malloc_Safely((void**)&sits_info.protein_water_frc, sizeof(VECTOR)*sits_info.atom_numbers);
//	Cuda_Malloc_Safely((void**)&sits_info.protein_water_energy, sizeof(float)*sits_info.protein_atom_numbers);
//	Cuda_Malloc_Safely((void**)&sits_info.atom_energy, sizeof(float)*sits_info.atom_numbers);
//
//	Cuda_Malloc_Safely((void**)&sits_info.sum_of_water_water_energy, sizeof(float));
//	Cuda_Malloc_Safely((void**)&sits_info.sum_of_protein_water_energy, sizeof(float));
//	Cuda_Malloc_Safely((void**)&sits_info.sum_of_protein_protein_energy, sizeof(float));
//
//
//	sits_info.k_numbers = amber_sits_info.ntemp;
//	Cuda_Malloc_Safely((void**)&sits_info.beta_k, sizeof(float)*sits_info.k_numbers);
//	Cuda_Malloc_Safely((void**)&sits_info.NkExpBetakU, sizeof(float)*sits_info.k_numbers);
//	Cuda_Malloc_Safely((void**)&sits_info.Nk, sizeof(float)*sits_info.k_numbers);
//	Cuda_Malloc_Safely((void**)&sits_info.sum_a, sizeof(float));
//	Cuda_Malloc_Safely((void**)&sits_info.sum_b, sizeof(float));
//	Cuda_Malloc_Safely((void**)&sits_info.factor, sizeof(float)*2);
//	Reset_List << <1, 2 >> >(2, sits_info.factor, 1.0);
//
//	float *beta_lin;
//	Malloc_Safely((void**)&beta_lin, sizeof(float)*sits_info.k_numbers);
//	float temp_slope = (amber_sits_info.temph - amber_sits_info.templ) / sits_info.k_numbers;
//	for (int i = 0; i < sits_info.k_numbers; i = i + 1)
//	{
//		beta_lin[i] = amber_sits_info.templ + temp_slope * i;
//		beta_lin[i] = 1. / (const_parameter.kB * beta_lin[i]);
//	}
//	//printf("%e %e\n",beta_lin[0],beta_lin[1]);
//	cudaMemcpy(sits_info.beta_k, beta_lin, sizeof(float)*sits_info.k_numbers, cudaMemcpyHostToDevice);
//
//
//	//SITS迭代部分初始化
//	sits_info.reset = 1;
//	sits_info.record_count = 0;
//	sits_info.record_interval = amber_sits_info.record_interval;
//	sits_info.update_interval = amber_sits_info.update_interval;
//
//	Cuda_Malloc_Safely((void**)&sits_info.ene_recorded, sizeof(float));
//	
//	Cuda_Malloc_Safely((void**)&sits_info.gf, sizeof(float) * sits_info.k_numbers);
//
//	Cuda_Malloc_Safely((void**)&sits_info.gfsum, sizeof(float));
//
//	Cuda_Malloc_Safely((void**)&sits_info.log_weight, sizeof(float)* sits_info.k_numbers);
//
//	Cuda_Malloc_Safely((void**)&sits_info.log_mk_inverse, sizeof(float)* sits_info.k_numbers);
//
//	Cuda_Malloc_Safely((void**)&sits_info.log_norm_old, sizeof(float)* sits_info.k_numbers);
//	Cuda_Malloc_Safely((void**)&sits_info.log_norm, sizeof(float)* sits_info.k_numbers);
//	if (amber_sits_info.fb_rest == 0)
//	{
//		for (int i = 0; i < sits_info.k_numbers; i = i + 1)
//		{
//			beta_lin[i] = -FLT_MAX;
//		}
//	}
//	else
//	{
//		FILE * norm_read_file;
//		Open_File_Safely(&norm_read_file, amber_sits_info.norm_Init_File_Name, "r");
//		for (int i = 0; i < sits_info.k_numbers; i = i + 1)
//		{
//			fscanf(norm_read_file, "%f", &beta_lin[i]);
//		}
//		fclose(norm_read_file);
//	}
//	cudaMemcpy(sits_info.log_norm_old, beta_lin, sizeof(float)* sits_info.k_numbers, cudaMemcpyHostToDevice);
//	cudaMemcpy(sits_info.log_norm, beta_lin, sizeof(float)* sits_info.k_numbers, cudaMemcpyHostToDevice);
//
//	Cuda_Malloc_Safely((void**)&sits_info.log_pk, sizeof(float)* sits_info.k_numbers);
//
//	Cuda_Malloc_Safely((void**)&sits_info.log_nk_inverse, sizeof(float)* sits_info.k_numbers);
//	for (int i = 0; i < sits_info.k_numbers; i = i + 1)
//	{
//		beta_lin[i] = 0.0;
//	}
//	cudaMemcpy(sits_info.log_nk_inverse, beta_lin, sizeof(float)* sits_info.k_numbers, cudaMemcpyHostToDevice);
//	
//	
//	Cuda_Malloc_Safely((void**)&sits_info.log_nk, sizeof(float)* sits_info.k_numbers);
//	if (amber_sits_info.fb_rest == 0)
//	{
//		for (int i = 0; i < sits_info.k_numbers; i = i + 1)
//		{
//			beta_lin[i] = amber_sits_info.fb_init;
//		}
//	}
//	else
//	{
//		FILE * fb_read_file;
//		Open_File_Safely(&fb_read_file,amber_sits_info.fb_Init_File_Name, "r");
//		for (int i = 0; i < sits_info.k_numbers; i = i + 1)
//		{
//			fscanf(fb_read_file, "%f", &beta_lin[i]);
//			beta_lin[i] = logf(beta_lin[i]);
//		}
//		fclose(fb_read_file);
//	}
//	cudaMemcpy(sits_info.log_nk, beta_lin, sizeof(float)* sits_info.k_numbers, cudaMemcpyHostToDevice);
//
//	//Nk的malloc已经在上面完成了
//	for (int i = 0; i < sits_info.k_numbers; i = i + 1)
//	{
//		beta_lin[i] = expf(beta_lin[i]);
//	}
//	cudaMemcpy(sits_info.Nk, beta_lin, sizeof(float)* sits_info.k_numbers, cudaMemcpyHostToDevice);
//
//
//	free(beta_lin);
//
//	
//
//	FILE *nk_record_file, *nk_rest_file, *norm_record_file, *norm_rest_file;
//	Open_File_Safely(&nk_record_file, amber_sits_info.fb_Record_File_Name, "wb");
//	Open_File_Safely(&nk_rest_file, amber_sits_info.fb_Rest_File_Name, "w");
//	Open_File_Safely(&norm_record_file, amber_sits_info.norm_Record_File_Name, "wb");
//	Open_File_Safely(&norm_rest_file, amber_sits_info.norm_Rest_File_Name, "w");
//
//	float *nk_recorded_cpu, *log_norm_recorded_cpu;
//	Malloc_Safely((void**)&nk_recorded_cpu, sizeof(float)*sits_info.k_numbers);
//	Malloc_Safely((void**)&log_norm_recorded_cpu, sizeof(float)*sits_info.k_numbers);
//
//
//	
//	/****************************************************************************************************/
//
//	Reset_List << <1, 1 >> >(1, nb_info.is_need_refresh_neighbor_list, 1);
//	SITS_Refresh_Neighbor_List << <1, 1 >> >
//		(nb_info.is_need_refresh_neighbor_list, 32,
//		md_info.atom_numbers, sits_info.protein_atom_numbers,
//		md_info.crd, md_info.old_crd, md_info.uint_crd,
//		0.5*md_info.crd_to_uint_crd_cof, md_info.uint_dr_to_dr_cof,
//		grid_info.atom_in_grid_serial,
//		2., md_info.box_length,
//		grid_info, grid_info.gpointer,
//		grid_info.bucket, grid_info.atom_numbers_in_grid_bucket,
//		nb_info.d_nl, nb_info.d_atomnumbers_in_nl,
//		sits_info.protein_water_neighbor, sits_info.protein_water_numbers,
//		nb_info.excluded_list_start, nb_info.excluded_list, nb_info.excluded_numbers);
//	
//
//
//	FILE *out, *ncFile;
//	Open_File_Safely(&out, amber_info.Out_File_Name, "w");
//	Open_File_Safely(&ncFile, amber_info.Traj_File_Name, "wb");
//
//
//	//主要循环
////	int refresh_count = 0;
//	dim3 thread_LJ(8, 32);
//	printf("start\n");
//	LARGE_INTEGER time_record;
//	LARGE_INTEGER m_nFreq;
//	QueryPerformanceFrequency(&m_nFreq);
//	QueryPerformanceCounter(&time_record);
//	int step_limit = amber_info.nstlim;
//	for (int steps = 1; steps <= step_limit; steps = steps + 1)
//	{
//		Crd_To_Uint_Crd << <ceilf((float)md_info.atom_numbers / 128), 128 >> >
//			(md_info.atom_numbers, 0.5*md_info.crd_to_uint_crd_cof, md_info.crd, md_info.uint_crd);
//		Copy_Crd_To_New_Crd << <ceilf((float)md_info.atom_numbers / 128), 128 >> >
//			(md_info.atom_numbers, md_info.uint_crd, md_info.uint_crd_with_LJ);
//
//		Reset_List << <ceilf((float)md_info.atom_numbers / 32), 32 >> >
//			(md_info.atom_numbers, sits_info.atom_energy, 0.);
//		Reset_List << <ceilf((float)sits_info.protein_atom_numbers / 32), 32 >> >
//			(sits_info.protein_atom_numbers, sits_info.protein_water_energy, 0.);
//		Reset_List << <ceilf((float)3.*md_info.atom_numbers / 32), 32 >> >
//			(3*md_info.atom_numbers, (float*)sits_info.protein_water_frc, 0.);
//
//		SITS_LJ_Force_With_Direct_CF_Energy << <ceilf((float)md_info.atom_numbers / 8), thread_LJ >> >
//			(md_info.atom_numbers, nb_info.d_nl, nb_info.d_atomnumbers_in_nl,
//			md_info.uint_crd_with_LJ, md_info.uint_dr_to_dr_cof,
//			LJ_info.d_LJ_A, LJ_info.d_LJ_B, nb_info.cutoff,
//			md_info.frc, pme_info.beta, 2. / sqrtf(const_parameter.pi),
//			sits_info.atom_energy);
//
//		SITS_Bond_Force_Energy << <ceilf((float)bond_info.bond_numbers / 128), 128 >> >
//			(bond_info.bond_numbers, md_info.uint_crd, md_info.uint_dr_to_dr_cof,
//			bond_info.atom_a, bond_info.atom_b, bond_info.k, bond_info.r0, md_info.frc, sits_info.atom_energy);
//
//
//		SITS_Angle_Force_Energy << <ceilf((float)angle_info.angle_numbers / 128), 128 >> >
//			(angle_info.angle_numbers,
//			md_info.uint_crd, md_info.uint_dr_to_dr_cof,
//			angle_info.atom_a, angle_info.atom_b, angle_info.atom_c,
//			angle_info.k, angle_info.theta0,
//			md_info.frc, sits_info.atom_energy);
//
//		SITS_Dihedral_Force_Energy << <ceilf((float)dihedral_info.dihedral_numbers / 128), 128 >> >(dihedral_info.dihedral_numbers, md_info.uint_crd, md_info.uint_dr_to_dr_cof,
//			dihedral_info.atom_a, dihedral_info.atom_b, dihedral_info.atom_c, dihedral_info.atom_d, dihedral_info.ipn, dihedral_info.pk, dihedral_info.gamc, dihedral_info.gams, dihedral_info.pn, md_info.frc,
//			sits_info.atom_energy);
//
//		SITS_Dihedral_14_LJ_Force_With_Direct_CF_Energy << <ceilf((float)dihedral_info.dihedral_14_numbers / 128), 128 >> >(dihedral_info.dihedral_14_numbers, md_info.uint_crd_with_LJ, md_info.uint_dr_to_dr_cof,
//			dihedral_info.atom_a_14, dihedral_info.atom_b_14, dihedral_info.lj_scale_factor, dihedral_info.cf_scale_factor, LJ_info.d_LJ_A, LJ_info.d_LJ_B, md_info.frc,
//			sits_info.atom_energy);
//
//		SITS_LJ_Force_With_Direct_CF_Energy << <ceilf((float)sits_info.protein_atom_numbers / 8), thread_LJ >> >
//			(sits_info.protein_atom_numbers, sits_info.protein_water_neighbor, sits_info.protein_water_numbers,
//			md_info.uint_crd_with_LJ, md_info.uint_dr_to_dr_cof,
//			LJ_info.d_LJ_A, LJ_info.d_LJ_B, nb_info.cutoff,
//			sits_info.protein_water_frc, pme_info.beta, 2. / sqrtf(const_parameter.pi),
//			sits_info.protein_water_energy);
//
//
//
//		Sum_Of_List << <1, 1024 >> >(0, sits_info.protein_atom_numbers, sits_info.atom_energy, sits_info.sum_of_protein_protein_energy);
//		Sum_Of_List << <1, 1024 >> >(sits_info.protein_atom_numbers, sits_info.atom_numbers, sits_info.atom_energy, sits_info.sum_of_water_water_energy);
//		Sum_Of_List << <1, 1024 >> >(sits_info.protein_atom_numbers, sits_info.protein_water_energy, sits_info.sum_of_protein_water_energy);
//
//		
//
//		if (amber_sits_info.fb_const == 0 && steps % sits_info.record_interval == 0)
//		{
//			SITS_Record_Ene << <1, 1 >> >(sits_info.ene_recorded, sits_info.sum_of_protein_water_energy, sits_info.sum_of_protein_protein_energy, amber_sits_info.pe_a, amber_sits_info.pe_b);
//			
//			SITS_Update_gf << <ceilf((float)sits_info.k_numbers / 32.), 32 >> >(sits_info.k_numbers, sits_info.gf,
//				sits_info.ene_recorded, sits_info.log_nk, sits_info.beta_k);
//
//			SITS_Update_gfsum << <1, 1 >> >(sits_info.k_numbers, sits_info.gfsum, sits_info.gf);
//			
//			
//			SITS_Update_log_pk << <ceilf((float)sits_info.k_numbers / 32.), 32 >> >(sits_info.k_numbers, sits_info.log_pk,
//				sits_info.gf, sits_info.gfsum, sits_info.reset);
//			
//			sits_info.reset = 0;
//			sits_info.record_count++;
//
//			if (sits_info.record_count % sits_info.update_interval == 0)
//			{
//				SITS_Update_log_mk_inverse << <ceilf((float)sits_info.k_numbers / 32.), 32 >> >(sits_info.k_numbers,
//					sits_info.log_weight, sits_info.log_mk_inverse, sits_info.log_norm_old,
//					sits_info.log_norm, sits_info.log_pk, sits_info.log_nk);
//
//				SITS_Update_log_nk_inverse << <1, 1 >> >(sits_info.k_numbers,
//					sits_info.log_nk_inverse, sits_info.log_mk_inverse);
//
//				SITS_Update_nk << <ceilf((float)sits_info.k_numbers / 32.), 32 >> >(sits_info.k_numbers,
//					sits_info.log_nk, sits_info.Nk, sits_info.log_nk_inverse);
//
//				
//				sits_info.record_count = 0;
//				sits_info.reset = 1;
//				
//
//				cudaMemcpy(nk_recorded_cpu, sits_info.Nk, sizeof(float)*sits_info.k_numbers, cudaMemcpyDeviceToHost);
//				fwrite(nk_recorded_cpu, sizeof(float), sits_info.k_numbers, nk_record_file);
//				cudaMemcpy(log_norm_recorded_cpu, sits_info.log_norm, sizeof(float)*sits_info.k_numbers, cudaMemcpyDeviceToHost);
//				fwrite(log_norm_recorded_cpu, sizeof(float), sits_info.k_numbers, norm_record_file);
//
//			}
//		}
//
//		
//
//		SITS_Enhanced_Force << <1, 1 >> >
//			(md_info.atom_numbers, sits_info.protein_atom_numbers,
//			md_info.frc, sits_info.protein_water_frc,
//			sits_info.sum_of_protein_protein_energy, sits_info.sum_of_protein_water_energy,
//			sits_info.k_numbers, sits_info.NkExpBetakU,
//			sits_info.beta_k, sits_info.Nk,
//			sits_info.sum_a, sits_info.sum_b, sits_info.factor,
//			1. / (liujian_info.target_temperature*const_parameter.kB), amber_sits_info.pe_a,amber_sits_info.pe_b,amber_sits_info.fb_bias);
//
//		
//
//		//Add_List << <1, 1024 >> >(3*md_info.atom_numbers,(float*) md_info.frc, (float*)sits_info.protein_water_frc);
//		pme_info.PME_Reciprocal_Force(md_info.uint_crd, md_info.d_charge, md_info.frc);
//
//		pme_info.PME_Excluded_Force(md_info.uint_crd, md_info.uint_dr_to_dr_cof,
//			md_info.d_charge, pme_info.beta, 2. / sqrtf(const_parameter.pi),
//			nb_info.excluded_list_start, nb_info.excluded_list, nb_info.excluded_numbers,
//			md_info.frc);
//
//
//		/*cudaMemcpy(md_info.force, md_info.frc, sizeof(VECTOR)*md_info.atom_numbers, cudaMemcpyDeviceToHost);
//		FILE *frc_record = fopen("20200722_Mannas_sys\\frc.txt", "w");
//		for (int i = 0; i < md_info.atom_numbers; i = i + 1)
//		{
//		fprintf(frc_record, "%f %f %f\n", md_info.force[i].x, md_info.force[i].y, md_info.force[i].z);
//		}
//		fclose(frc_record);
//		printf("done");
//		getchar(); */
//
//		if (steps % amber_info.ntwx == 0)
//		{
//			float s = 0., lj_ene = 0., bond_ene = 0., angle_ene = 0.;
//
//			MD_Residue_Ek << <ceilf((float)md_info.residue_numbers / 32.), 32 >> >
//				(md_info.residue_numbers, md_info.res_start, md_info.res_end, md_info.res_ek_energy,
//				md_info.vel, md_info.d_mass);
//			Sum_Of_List << <1, 1024 >> >(md_info.residue_numbers, md_info.res_ek_energy, md_info.sigma_of_res_ek);
//
//
//
//			LJ_info.LJ_Energy(md_info.uint_crd_with_LJ, md_info.uint_dr_to_dr_cof,
//				nb_info.d_nl, nb_info.d_atomnumbers_in_nl, nb_info.cutoff_square);
//
//			Bond_Energy << <ceilf((float)bond_info.bond_numbers / 32), 32 >> >
//				(bond_info.bond_numbers, md_info.uint_crd, md_info.uint_dr_to_dr_cof,
//				bond_info.atom_a, bond_info.atom_b, bond_info.k, bond_info.r0, bond_info.bond_ene);
//			Sum_Of_List << <1, 1024 >> >(bond_info.bond_numbers, bond_info.bond_ene, bond_info.sigme_of_bond_ene);
//
//			Angle_Energy << <ceilf((float)angle_info.angle_numbers / 32), 32 >> >
//				(angle_info.angle_numbers,
//				md_info.uint_crd, md_info.uint_dr_to_dr_cof,
//				angle_info.atom_a, angle_info.atom_b, angle_info.atom_c,
//				angle_info.k, angle_info.theta0,
//				angle_info.energy);
//			Sum_Of_List << <1, 1024 >> >(angle_info.angle_numbers, angle_info.energy, angle_info.sigma_energy);
//
//			float dihedral_ene, dihedral_14_LJ_ene, dihedral_14_CF_ene;
//			Dihedral_Energy << <ceilf((float)dihedral_info.dihedral_numbers / 32), 32 >> >(dihedral_info.dihedral_numbers, md_info.uint_crd, md_info.uint_dr_to_dr_cof,
//				dihedral_info.atom_a, dihedral_info.atom_b, dihedral_info.atom_c, dihedral_info.atom_d, dihedral_info.ipn, dihedral_info.pk, dihedral_info.gamc, dihedral_info.gams, dihedral_info.pn, dihedral_info.dihedral_energy);
//			Sum_Of_List << <1, 1024 >> >
//				(dihedral_info.dihedral_numbers, dihedral_info.dihedral_energy, dihedral_info.sigma_energy);
//
//			Dihedral_14_LJ_Energy << <ceilf((float)dihedral_info.dihedral_14_numbers / 32), 32 >> >
//				(dihedral_info.dihedral_14_numbers, md_info.uint_crd_with_LJ, md_info.uint_dr_to_dr_cof,
//				dihedral_info.atom_a_14, dihedral_info.atom_b_14, dihedral_info.lj_scale_factor,
//				LJ_info.d_LJ_A, LJ_info.d_LJ_B, dihedral_info.dihedral_14_energy);
//			Sum_Of_List << <1, 1024 >> >
//				(dihedral_info.dihedral_14_numbers, dihedral_info.dihedral_14_energy, dihedral_info.sigma_14_LJ_energy);
//
//			Dihedral_14_CF_Energy << <ceilf((float)dihedral_info.dihedral_14_numbers / 32), 32 >> >
//				(dihedral_info.dihedral_14_numbers, md_info.uint_crd_with_LJ, md_info.uint_dr_to_dr_cof,
//				dihedral_info.atom_a_14, dihedral_info.atom_b_14, dihedral_info.cf_scale_factor,
//				dihedral_info.dihedral_14_energy);
//			Sum_Of_List << <1, 1024 >> >
//				(dihedral_info.dihedral_14_numbers, dihedral_info.dihedral_14_energy, dihedral_info.sigma_14_CF_energy);
//
//
//			pme_info.PME_Energy(md_info.uint_crd, md_info.d_charge,
//				nb_info.d_nl, nb_info.d_atomnumbers_in_nl, md_info.uint_dr_to_dr_cof,
//				nb_info.excluded_list_start, nb_info.excluded_list, nb_info.excluded_numbers);
//
//			
//			cudaMemcpy(&s, md_info.sigma_of_res_ek, sizeof(float), cudaMemcpyDeviceToHost);
//			cudaMemcpy(&lj_ene, LJ_info.d_LJ_energy_sum, sizeof(float), cudaMemcpyDeviceToHost);
//			cudaMemcpy(&bond_ene, bond_info.sigme_of_bond_ene, sizeof(float), cudaMemcpyDeviceToHost);
//			cudaMemcpy(&angle_ene, angle_info.sigma_energy, sizeof(float), cudaMemcpyDeviceToHost);
//			cudaMemcpy(&dihedral_ene, dihedral_info.sigma_energy, sizeof(float), cudaMemcpyDeviceToHost);
//			cudaMemcpy(&dihedral_14_CF_ene, dihedral_info.sigma_14_CF_energy, sizeof(float), cudaMemcpyDeviceToHost);
//			cudaMemcpy(&dihedral_14_LJ_ene, dihedral_info.sigma_14_LJ_energy, sizeof(float), cudaMemcpyDeviceToHost);
//			float tempFb;
//			cudaMemcpy(&tempFb, sits_info.factor, sizeof(float), cudaMemcpyDeviceToHost);
//			printf("_steps_ _TEMP_ _LJ_ENE_ _BOND_ENE_ _ANGLE_ENE_ _DIHEDRAL_ENE_ _14LJ_ENE_ _14CF_ENE_ _CF_PME_ENE_\n");
//			printf("%7d %6.2f %6.2f   %6.2f     %6.2f ",
//				steps, 2. / 3. / const_parameter.kB / md_info.residue_numbers*s, lj_ene, bond_ene, angle_ene);
//			printf("       %6.2f     %6.2f  %6.2f        %6.2f\n", dihedral_ene, dihedral_14_LJ_ene, dihedral_14_CF_ene, pme_info.ee_ene);
//			fprintf(out,"%7d %6.2f %6.2f   %6.2f     %6.2f ",
//				steps, 2. / 3. / const_parameter.kB / md_info.residue_numbers*s, lj_ene, bond_ene, angle_ene);
//			fprintf(out,"       %6.2f     %6.2f  %6.2f        %6.2f\n", dihedral_ene, dihedral_14_LJ_ene, dihedral_14_CF_ene, pme_info.ee_ene);
//			printf("__fb__\n%6.2f\n", tempFb+amber_sits_info.fb_bias);
//			printf("\n");
//
//			
//			float f,ff,fff;
//			cudaMemcpy(&f, sits_info.sum_of_water_water_energy, sizeof(float), cudaMemcpyDeviceToHost);
//			cudaMemcpy(&ff, sits_info.sum_of_protein_water_energy, sizeof(float), cudaMemcpyDeviceToHost);
//			cudaMemcpy(&fff, sits_info.sum_of_protein_protein_energy, sizeof(float), cudaMemcpyDeviceToHost);
//			printf("sits energy test %f %f %f %f %f %f\n", f + ff + fff, f, ff, fff, (float)(fff + 0.5*ff) ,
//				angle_ene + bond_ene + dihedral_ene + dihedral_14_LJ_ene + dihedral_14_CF_ene + pme_info.ee_ene + lj_ene);
//
//			fprintf(out, "%d %f %f %f %f\n", steps, fff, ff, f, tempFb + amber_sits_info.fb_bias);
//
//			Vector_Translation << <ceilf((float)md_info.atom_numbers / 32), 32 >> >(md_info.atom_numbers, md_info.crd, trans_vec_minus);
//			cudaMemcpy(md_info.coordinate, md_info.crd, sizeof(VECTOR)*md_info.atom_numbers, cudaMemcpyDeviceToHost);
//			Vector_Translation << <ceilf((float)md_info.atom_numbers / 32), 32 >> >(md_info.atom_numbers, md_info.crd, trans_vec);
//			fwrite(&md_info.coordinate[0].x, sizeof(VECTOR), md_info.atom_numbers, ncFile);
//			cudaMemcpy(md_info.velocity, md_info.vel, sizeof(VECTOR)*md_info.atom_numbers, cudaMemcpyDeviceToHost);
//			//Export_Information_To_Rst7(amber_info.Rst_File_Name, amber_info.md_name, md_info.simulation_start_time + amber_info.dt * steps / 1000.0, md_info.atom_numbers, md_info.coordinate, md_info.velocity, md_info.box_length);
//		}
//
//		/*MD_Iteration_Leap_Frog << <ceilf((float)md_info.atom_numbers / 128), 128 >> >
//			(md_info.atom_numbers, md_info.vel, md_info.crd, md_info.frc, md_info.acc, md_info.d_mass_inverse, md_info.dt);*/
//		Rand_Normal << <ceilf((float)liujian_info.float4_numbers / 32.), 32 >> >
//			(liujian_info.float4_numbers, liujian_info.rand_state, (float4 *)liujian_info.random_force);
//		MD_Iteration_Leap_Frog_With_LiuJian
//			<< <ceilf((float)md_info.atom_numbers / 32), 32 >> >
//			(md_info.atom_numbers, 0.5*md_info.dt, md_info.dt, liujian_info.exp_gamma, md_info.d_mass_inverse,
//			liujian_info.sqrt_mass,md_info.vel, md_info.crd, md_info.frc, md_info.acc, liujian_info.random_force);
//
//
//		Is_need_refresh_neighbor_list_cuda << <ceilf((float)md_info.atom_numbers / 128), 128 >> >
//			(md_info.atom_numbers, md_info.crd, md_info.old_crd, 1., nb_info.is_need_refresh_neighbor_list);
//		/*Refresh_Neighbor_List << <1, 1 >> >
//			(nb_info.is_need_refresh_neighbor_list, 32,
//			md_info.atom_numbers, md_info.crd, md_info.old_crd, md_info.uint_crd,
//			0.5*md_info.crd_to_uint_crd_cof, md_info.uint_dr_to_dr_cof,
//			grid_info.atom_in_grid_serial,
//			2., md_info.box_length,
//			grid_info, grid_info.gpointer,
//			grid_info.bucket, grid_info.atom_numbers_in_grid_bucket,
//			nb_info.d_nl, nb_info.d_atomnumbers_in_nl, nb_info.excluded_list_start, nb_info.excluded_list, nb_info.excluded_numbers);*/
//
//		SITS_Refresh_Neighbor_List << <1, 1 >> >
//			(nb_info.is_need_refresh_neighbor_list, 32,
//			md_info.atom_numbers, sits_info.protein_atom_numbers,
//			md_info.crd, md_info.old_crd, md_info.uint_crd,
//			0.5*md_info.crd_to_uint_crd_cof, md_info.uint_dr_to_dr_cof,
//			grid_info.atom_in_grid_serial,
//			2., md_info.box_length,
//			grid_info, grid_info.gpointer,
//			grid_info.bucket, grid_info.atom_numbers_in_grid_bucket,
//			nb_info.d_nl, nb_info.d_atomnumbers_in_nl,
//			sits_info.protein_water_neighbor, sits_info.protein_water_numbers,
//			nb_info.excluded_list_start, nb_info.excluded_list, nb_info.excluded_numbers);
//
//
//	}
//
//	cudaMemcpy(nk_recorded_cpu, sits_info.Nk, sizeof(float)*sits_info.k_numbers, cudaMemcpyDeviceToHost);
//	cudaMemcpy(log_norm_recorded_cpu, sits_info.log_norm, sizeof(float)*sits_info.k_numbers, cudaMemcpyDeviceToHost);
//	for (int i = 0; i < sits_info.k_numbers; i++)
//	{
//		fprintf(nk_rest_file,"%e\n", nk_recorded_cpu[i]);
//		fprintf(norm_rest_file, "%e\n", log_norm_recorded_cpu[i]);
//	}
//	Vector_Translation << <ceilf((float)md_info.atom_numbers / 32), 32 >> >(md_info.atom_numbers, md_info.crd, trans_vec_minus);
//	cudaMemcpy(md_info.coordinate, md_info.crd, sizeof(VECTOR)*md_info.atom_numbers, cudaMemcpyDeviceToHost);
//	cudaMemcpy(md_info.velocity, md_info.vel, sizeof(VECTOR)*md_info.atom_numbers, cudaMemcpyDeviceToHost);
//	Export_Information_To_Rst7(amber_info.Rst_File_Name, amber_info.md_name, amber_info.dt * step_limit, md_info.atom_numbers, md_info.coordinate, md_info.velocity, md_info.box_length);
//	LARGE_INTEGER time_record2;
//	QueryPerformanceCounter(&time_record2);
//	float total_time = (float)(time_record2.QuadPart - time_record.QuadPart) / m_nFreq.QuadPart;
//	printf("time %f\nspeed %f ns/day\n", total_time, step_limit * amber_info.dt / 1000000 / (total_time / 86400.0));	
//	fclose(out);
//	fclose(nk_record_file); 
//	fclose(nk_rest_file);
//	fclose(ncFile);
//	fclose(norm_record_file);
//	fclose(norm_record_file);
//	return 0;
//}
//
//int main(int argc, char * argv[])
//{
//	AMBER_SITS_INFORMATION amber_sits_info;
//	AMBER_INFORMATION amber_info;
//	strcpy(amber_info.Rst_File_Name, "restrt");
//	strcpy(amber_info.Out_File_Name, "mdout");
//	strcpy(amber_info.Traj_File_Name, "mdcrd");
//	strcpy(amber_info.Parm_File_Name, "parm7");
//	strcpy(amber_info.Crd_File_Name, "rst7");
//	strcpy(amber_info.In_File_Name, "mdin");
//	//strcpy(amber_info.Parm_File_Name, "sits_test//ala_wat.parm7");
//	//strcpy(amber_info.Crd_File_Name, "sits_test//3.rst7");
//	//strcpy(amber_info.In_File_Name, "sits_test//test.in");
//	for (int i = 0; i < argc; i++)
//	{
//		if (strcmp(argv[i], "-i") == 0)
//		{
//			i++;
//
//			if (strlen(argv[i]) > 255)
//			{
//				printf("%s is too long.\nfile name can not be longer than 255 bits.", argv[i]);
//				getchar();
//				exit(1);
//			}
//
//			strcpy(amber_info.In_File_Name, argv[i]);
//		}
//		else if (strcmp(argv[i], "-p") == 0)
//		{
//			i++;
//
//			if (strlen(argv[i]) > 255)
//			{
//				printf("%s is too long.\nfile name can not be longer than 255 bits.", argv[i]);
//				getchar();
//				exit(1);
//			}
//
//			strcpy(amber_info.Parm_File_Name, argv[i]);
//		}
//		else if (strcmp(argv[i], "-c") == 0)
//		{
//			i++;
//
//			if (strlen(argv[i]) > 255)
//			{
//				printf("%s is too long.\nfile name can not be longer than 255 bits.", argv[i]);
//				getchar();
//				exit(1);
//			}
//
//			strcpy(amber_info.Crd_File_Name, argv[i]);
//		}
//		else if (strcmp(argv[i], "-r") == 0)
//		{
//			i++;
//
//			if (strlen(argv[i]) > 255)
//			{
//				printf("%s is too long.\nfile name can not be longer than 255 bits.", argv[i]);
//				getchar();
//				exit(1);
//			}
//
//			strcpy(amber_info.Rst_File_Name, argv[i]);
//		}
//		else if (strcmp(argv[i], "-x") == 0)
//		{
//			i++;
//
//			if (strlen(argv[i]) > 255)
//			{
//				printf("%s is too long.\nfile name can not be longer than 255 bits.", argv[i]);
//				getchar();
//				exit(1);
//			}
//
//			strcpy(amber_info.Traj_File_Name, argv[i]);
//		}
//		else if (strcmp(argv[i], "-o") == 0)
//		{
//			i++;
//
//			if (strlen(argv[i]) > 255)
//			{
//				printf("%s is too long.\nfile name can not be longer than 255 bits.", argv[i]);
//				getchar();
//				exit(1);
//			}
//
//			strcpy(amber_info.Out_File_Name, argv[i]);
//		}
//	}
//	amber_sits_info.amber_info = amber_info;
//	printf("FILE LIST\n  Input Control: %s\n  Start Coordinates: %s\n  Topology: %s\n", amber_sits_info.amber_info.In_File_Name, amber_info.Crd_File_Name, amber_info.Parm_File_Name);
//	printf("  Restart Coordinates: %s\n  Coordinates Trajectory: %s\n  Energy Out: %s\n", amber_info.Rst_File_Name, amber_info.Traj_File_Name, amber_info.Out_File_Name);
//
//	Import_MD_Control_Information_From_AMBERFILE(&amber_sits_info);
//
//	printf("SITS FILE LIST\n  norm_Record_File_Name: %s\n  norm_Rest_File_Name: %s\n  norm_Init_File_Name: %s\n", amber_sits_info.norm_Record_File_Name, amber_sits_info.norm_Rest_File_Name, amber_sits_info.norm_Init_File_Name);
//	printf("  fb_Record_File_Name: %s\n  fb_Rest_File_Name: %s\n  fb_Init_File_Name: %s\n", amber_sits_info.fb_Record_File_Name, amber_sits_info.fb_Rest_File_Name, amber_sits_info.fb_Init_File_Name);
//
//	int exit_code;
//	if (amber_sits_info.amber_info.mode == 1)
//	{
//		exit_code = main_SITS(amber_sits_info);
//	}
//
//
//	if (exit_code == 0)
//	{
//		printf("md end!\n");
//		//getchar();
//		return 0;
//	}
//}