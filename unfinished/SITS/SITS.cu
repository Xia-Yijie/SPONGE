#include "SITS.cuh"

static __global__ void SITS_For_Enhanced_Force_Protein
(const int protein_numbers, VECTOR *md_frc, const VECTOR *pw_frc, const float fc_ball, const float factor)
{
	for (int i = threadIdx.x; i < protein_numbers; i = i + blockDim.x)
	{
		md_frc[i].x = fc_ball*(md_frc[i].x) + factor*pw_frc[i].x;
		md_frc[i].y = fc_ball*(md_frc[i].y) + factor*pw_frc[i].y;
		md_frc[i].z = fc_ball*(md_frc[i].z) + factor*pw_frc[i].z;
	}
}
static __global__ void SITS_For_Enhanced_Force_Water
(const int protein_numbers, const int atom_numbers, VECTOR *md_frc, const VECTOR *pw_frc, const float factor)
{

	for (int i = threadIdx.x + protein_numbers; i < atom_numbers; i = i + blockDim.x)
	{
		md_frc[i].x = md_frc[i].x + factor*pw_frc[i].x;
		md_frc[i].y = md_frc[i].y + factor*pw_frc[i].y;
		md_frc[i].z = md_frc[i].z + factor*pw_frc[i].z;
	}
}

static __device__ __host__ float log_add_log(float a, float b)
{
	return fmaxf(a, b) + logf(1.0 + expf(-fabsf(a - b)));
}

static __global__ void SITS_Record_Ene(float *ene_record, const float *pw_ene, const float *pp_ene,
	const float pe_a, const float pe_b, const float pwwp_factor)
{
	float temp = *pw_ene * pwwp_factor + *pp_ene;
	temp = pe_a * temp + pe_b;
	
	*ene_record = temp;
	//printf("DEBUG ene_record: %f\n", ene_record[0]);
}

static __global__ void SITS_Update_gf(const int kn, float *gf, 
	const float *ene_record, const float *log_nk, const float *beta_k)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < kn)
	{
		gf[i] = -beta_k[i] * ene_record[0] + log_nk[i];
		//printf("DEBUG gf: %d %f\n", i, gf[i]);
	}
}

static __global__ void SITS_Update_gfsum(const int kn, float *gfsum, const float *gf)
{
	if (threadIdx.x == 0)
	{
		gfsum[0] = -FLT_MAX;
	}
	for (int i = 0; i < kn; i = i + 1)
	{
		gfsum[0] = log_add_log(gfsum[0], gf[i]);
		//printf("DEBUG gfsum: %d %f %f\n", i, gfsum[0], gf[i]);
	}
}

static __global__ void SITS_Update_log_pk(const int kn, float *log_pk, 
	const float *gf, const float *gfsum, const int reset)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < kn)
	{
		float gfi = gf[i];
		log_pk[i] = ((float)reset) * gfi + ((float)(1 - reset)) * log_add_log(log_pk[i], gfi - gfsum[0]);
		//printf("DEBUG log_pk: %d %f %f\n", i, log_pk[i], gfsum[0]);
	}
}


static __global__ void SITS_Update_log_mk_inverse(const int kn, 
	float *log_weight, float *log_mk_inverse, float *log_norm_old, 
	float *log_norm, const float *log_pk, const float *log_nk)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < kn - 1)
	{
		log_weight[i] = (log_pk[i] + log_pk[i + 1]) * 0.5;
		//printf("DEBUG log_weight: %d %f %f\n", i, log_pk[i], log_pk[i + 1]);
		log_mk_inverse[i] = log_nk[i] - log_nk[i + 1];
		log_norm_old[i] = log_norm[i];
		log_norm[i] = log_add_log(log_norm[i], log_weight[i]);
		log_mk_inverse[i] = log_add_log(log_mk_inverse[i] + log_norm_old[i] - log_norm[i], log_pk[i + 1] - log_pk[i] + log_mk_inverse[i] + log_weight[i] - log_norm[i]);
		//printf("DEBUG log_norm: %d %f %f\n", i, log_norm[i], log_weight[i]);
	}
}

static __global__ void SITS_Update_log_nk_inverse(const int kn,
	float *log_nk_inverse, 	const float *log_mk_inverse)
{
	for (int i = 0; i < kn - 1; i++)
	{
		log_nk_inverse[i + 1] = log_nk_inverse[i] + log_mk_inverse[i];
	}
}

static __global__ void SITS_Update_nk(const int kn,
	float *log_nk, float *nk, const float *log_nk_inverse)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < kn )
	{
		log_nk[i] = -log_nk_inverse[i];
		nk[i] = exp(log_nk[i]);
	}
}

__global__ void SITS_For_Enhanced_Force_Calculate_NkExpBetakU_1
(const int k_numbers, const float *beta_k, const float *nk,float *nkexpbetaku,
const float ene)
{
	float lin = beta_k[k_numbers-1];
	for (int i = threadIdx.x; i < k_numbers; i = i + blockDim.x)
	{
		nkexpbetaku[i] = nk[i] * expf(-(beta_k[i] - lin) * ene);
		//printf("%f %f\n", beta_k[i], nkexpbetaku[i]);
	}
}

__global__ void SITS_For_Enhanced_Force_Calculate_NkExpBetakU_2
(const int k_numbers, const float *beta_k, const float *nk, float *nkexpbetaku,
const float ene)
{
	float lin = beta_k[0];
	for (int i = threadIdx.x; i < k_numbers; i = i + blockDim.x)
	{
		nkexpbetaku[i] = nk[i] * expf(-(beta_k[i] - lin) * ene);
		//printf("%f %f\n", beta_k[i], nkexpbetaku[i]);
	}
}

__global__ void SITS_For_Enhanced_Force_Sum_Of_Above
(const int k_numbers, const float *nkexpbetaku, const float *beta_k, float *sum_of_above)
{
	if (threadIdx.x == 0)
	{
		sum_of_above[0] = 0.;
	}
	__syncthreads();
	float lin = 0.;
	for (int i = threadIdx.x; i < k_numbers; i = i + blockDim.x)
	{
		lin = lin + beta_k[i]*nkexpbetaku[i];
	}
	atomicAdd(sum_of_above, lin);
}

__global__ void SITS_For_Enhanced_Force_Sum_Of_NkExpBetakU
(const int k_numbers,const float *nkexpbetaku, float *sum_of_below)
{
	if (threadIdx.x == 0)
	{
		sum_of_below[0] = 0.;
	}
	__syncthreads();
	float lin = 0.;
	for (int i = threadIdx.x; i < k_numbers; i = i + blockDim.x)
	{
		lin = lin + nkexpbetaku[i];
		//printf("%f\n", nkexpbetaku[i]);
	}
	atomicAdd(sum_of_below, lin);
}

__global__ void SITS_Classical_Enhanced_Force
(const int atom_numbers, const int protein_atom_numbers, const float pwwp_factor,
VECTOR *md_frc, const VECTOR *pw_frc,
const float *pp_ene, const float *pw_ene,
const int k_numbers,float *nkexpbetaku,
const float *beta_k, const float *n_k,
float *sum_a,float *sum_b,float *factor,
const float beta_0, const float pe_a, const float pe_b, const float fb_bias)
{
	float ene = pp_ene[0] + pwwp_factor * pw_ene[0];
	ene = pe_a * ene + pe_b;
	if (ene > 0)
	{
		SITS_For_Enhanced_Force_Calculate_NkExpBetakU_1 << <1, 64 >> >
			(k_numbers, beta_k, n_k, nkexpbetaku, ene);
	}
	else
	{
		SITS_For_Enhanced_Force_Calculate_NkExpBetakU_2 << <1, 64 >> >
			(k_numbers, beta_k, n_k, nkexpbetaku, ene);
	}

	SITS_For_Enhanced_Force_Sum_Of_NkExpBetakU << <1, 128 >> >
			(k_numbers, nkexpbetaku, sum_b);

	SITS_For_Enhanced_Force_Sum_Of_Above << <1, 128 >> >
		(k_numbers, nkexpbetaku, beta_k, sum_a);


	
	factor[0] = sum_a[0] / sum_b[0] / beta_0 + fb_bias;
	//这段代码是避免fc_ball变化太大而造成体系崩溃的
	if (!isinf(factor[0]) && !isnan(factor[0]) && (factor[0] > 0.5 * factor[1]) && (factor[0] < 2 * factor[1]) )
	{
		factor[1] = factor[0];
	}
	else
	{
		factor[0] = factor[1];
	}
	float fc = factor[0];
	
//	printf("factor %e sum0 %e %e ene %f lfactor %e\n", fc, sum_a[0], sum_b[0], ene, factor[1]);
	__syncthreads();
	

	//line
	//fc = (ene - 20.) / 80./2. + 0.2;
	SITS_For_Enhanced_Force_Protein << <1, 128 >> >
		(protein_atom_numbers, md_frc, pw_frc, fc, pwwp_factor * fc + 1.0 - pwwp_factor);
	SITS_For_Enhanced_Force_Water << <1, 128 >> >
		(protein_atom_numbers, atom_numbers, md_frc, pw_frc, pwwp_factor * fc + 1.0 - pwwp_factor);
}



void SITS::SITS_Classical_Update_Info(int steps)
{
	if (!classical_info.constant_nk && steps % classical_info.record_interval == 0)
	{
		SITS_Record_Ene << <1, 1 >> >(classical_info.ene_recorded, 
			d_total_protein_water_atom_energy, d_total_pp_atom_energy, classical_info.energy_multiple, classical_info.energy_shift,
			info.pwwp_enhance_factor);
				
		SITS_Update_gf << <ceilf((float)classical_info.k_numbers / 32.), 32 >> >(classical_info.k_numbers, classical_info.gf,
			classical_info.ene_recorded, classical_info.log_nk, classical_info.beta_k);
	
		SITS_Update_gfsum << <1, 1 >> >(classical_info.k_numbers, classical_info.gfsum, classical_info.gf);
				
				
		SITS_Update_log_pk << <ceilf((float)classical_info.k_numbers / 32.), 32 >> >(classical_info.k_numbers, classical_info.log_pk,
			classical_info.gf, classical_info.gfsum, classical_info.reset);
				
		classical_info.reset = 0;
		classical_info.record_count++;
	
		if (classical_info.record_count % classical_info.update_interval == 0)
		{
			SITS_Update_log_mk_inverse << <ceilf((float)classical_info.k_numbers / 32.), 32 >> >(classical_info.k_numbers,
				classical_info.log_weight, classical_info.log_mk_inverse, classical_info.log_norm_old,
				classical_info.log_norm, classical_info.log_pk, classical_info.log_nk);
	
			SITS_Update_log_nk_inverse << <1, 1 >> >(classical_info.k_numbers,
				classical_info.log_nk_inverse, classical_info.log_mk_inverse);
	
			SITS_Update_nk << <ceilf((float)classical_info.k_numbers / 32.), 32 >> >(classical_info.k_numbers,
				classical_info.log_nk, classical_info.Nk, classical_info.log_nk_inverse);
	
					
			classical_info.record_count = 0;
			classical_info.reset = 1;
					
            if (!classical_info.constant_nk)
            {
                cudaMemcpy(classical_info.log_nk_recorded_cpu, classical_info.Nk, sizeof(float)*classical_info.k_numbers, cudaMemcpyDeviceToHost);
                fwrite(classical_info.log_nk_recorded_cpu, sizeof(float), classical_info.k_numbers, classical_info.nk_traj_file);
                cudaMemcpy(classical_info.log_norm_recorded_cpu, classical_info.log_norm, sizeof(float)*classical_info.k_numbers, cudaMemcpyDeviceToHost);
                fwrite(classical_info.log_norm_recorded_cpu, sizeof(float), classical_info.k_numbers, classical_info.norm_traj_file);
            }
		}
	}
}

void SITS::Initial_SITS(CONTROLLER *controller,int atom_numbers)
{
	printf("\nInitial Selective Integrated Tempering Sampling:\n");
	if (controller[0].Command_Exist("sits_mode"))
	{
        if (strcmp(controller[0].Command("sits_mode"), "only_select_force") == 0)
        {
            info.sits_mode = ONLY_SELECT_MODE_FORCE;
            printf("	Choosing No Enhancing Mode, Only select force.\n");
        }
        else if (strcmp(controller[0].Command("sits_mode"), "only_select_force_energy") == 0)
        {
            info.sits_mode = ONLY_SELECT_MODE_FORCE_ENERGY;
            printf("	Choosing No Enhancing Mode, Only select force and energy.\n");
        }
		else if (strcmp(controller[0].Command("sits_mode"), "simple") == 0)
		{
			info.sits_mode = SIMPLE_SITS_MODE;
			printf("	Choosing SimpleSITS Mode\n");
		}
		else if (strcmp(controller[0].Command("sits_mode"), "classical") == 0)
		{
			info.sits_mode = CLASSICAL_SITS_MODE;
			printf("	Choosing ClassicalSITS Mode\n");
		}
		else
		{
			printf("Error: Please Choose a Correct SITS Mode!\n");
			getchar();
			exit(1);
		}
	}
	else
	{
		printf("Error: Please Choose a SITS Mode!\n");
		getchar();
	}


	if (controller[0].Command_Exist("sits_atom_numbers"))
	{
		info.protein_atom_numbers = atoi(controller[0].Command("sits_atom_numbers"));
		printf("	sits_atom_numbers is %d\n", info.protein_atom_numbers);
	}
	else
	{
		printf("Error: Please Give a sits_atom_numbers!\n");
		getchar();
		exit(1);
	}
	info.atom_numbers = atom_numbers;

	//初始化一些数组
	info.max_neighbor_numbers = 800;
	if (controller[0].Command_Exist("max_neighbor_numbers"))
		info.max_neighbor_numbers = atoi(controller[0].Command("max_neighbor_numbers"));
	Cuda_Malloc_Safely((void**)&d_protein_water_atom_energy, sizeof(float)*info.atom_numbers);
	Cuda_Malloc_Safely((void**)&d_total_pp_atom_energy, sizeof(float)* 1);
	Cuda_Malloc_Safely((void**)&d_total_ww_atom_energy, sizeof(float)* 1);
	Cuda_Malloc_Safely((void**)&d_total_protein_water_atom_energy, sizeof(float)*1);
	Cuda_Malloc_Safely((void**)&protein_water_frc, sizeof(VECTOR)*info.atom_numbers);
	Reset_List << <ceilf((float)3.*info.atom_numbers / 128), 128 >> >(3 * info.atom_numbers, (float*)protein_water_frc, 0.);

	ATOM_GROUP *h_nl=NULL;//临时的近邻表
	Malloc_Safely((void**)&h_nl, sizeof(ATOM_GROUP)*info.atom_numbers);
	Cuda_Malloc_Safely((void**)&d_nl_ppww, sizeof(ATOM_GROUP)*info.atom_numbers);
	Cuda_Malloc_Safely((void**)&d_nl_pwwp, sizeof(ATOM_GROUP)*info.atom_numbers);
	for (int i = 0; i < info.atom_numbers; i = i + 1)
	{
		h_nl[i].atom_numbers = 0;
		Cuda_Malloc_Safely((void**)&h_nl[i].atom_serial, sizeof(int)* info.max_neighbor_numbers);
	}
	cudaMemcpy(d_nl_ppww, h_nl, sizeof(ATOM_GROUP)*info.atom_numbers, cudaMemcpyHostToDevice);
	
	for (int i = 0; i < info.atom_numbers; i = i + 1)
	{
		h_nl[i].atom_numbers = 0;
		Cuda_Malloc_Safely((void**)&h_nl[i].atom_serial, sizeof(int)* info.max_neighbor_numbers);
	}
	cudaMemcpy(d_nl_pwwp, h_nl, sizeof(ATOM_GROUP)*info.atom_numbers, cudaMemcpyHostToDevice);
	free(h_nl);
    

	//记录分能量的文件
	if (controller[0].Command_Exist("sits_energy_record"))
	{
        printf("	SITS Energy Record File: %s\n", controller[0].Command("sits_energy_record"));
		Open_File_Safely(&sits_ene_record_out, controller[0].Command("sits_energy_record"), "w");
	}
	else
	{
        printf("	SITS Energy Record File: SITS_Energy_Record.txt\n");
		Open_File_Safely(&sits_ene_record_out,"SITS_Energy_Record.txt", "w");
	}
    fprintf(sits_ene_record_out, "SITS: ______AA______ ______BB______ ______AB______ fc_ball\n");
	if (info.sits_mode % 2 == 1)
	{
		info.pwwp_enhance_factor = 0.5;
		if (controller[0].Command_Exist("pwwp_enhance_factor"))
		{
			info.pwwp_enhance_factor = atof(controller[0].Command("pwwp_enhance_factor"));
		}
		printf("	SITS Interaction Enhancing Factor: %.2f\n", info.pwwp_enhance_factor);
	}

	//simple SITS模式的额外初始化
	if (info.sits_mode == SIMPLE_SITS_MODE)
	{
		//初始化随机种子，用于fc随机运动
		srand(simple_info.random_seed);
		//初始化fc_pdf
		if (controller[0].Command_Exist("sits_fcball_pdf_grid_numbers"))
		{
			simple_info.grid_numbers = atoi(controller[0].Command("sits_fcball_pdf_grid_numbers"));
		}
		Malloc_Safely((void**)&simple_info.fc_pdf, sizeof(float)*simple_info.grid_numbers);
		for (int i = 0; i < simple_info.grid_numbers; i = i + 1)
		{
			simple_info.fc_pdf[i] = 0.01;//默认概率都为相同的非0值
		}
		if (controller[0].Command_Exist("sits_fcball_pdf"))
		{
			FILE *flin=NULL;//打开读入的临时文件指针
			Open_File_Safely(&flin,controller[0].Command("sits_fcball_pdf"),"r");
			for (int i = 0; i < simple_info.grid_numbers; i = i + 1)
			{
				fscanf(flin,"%f\n",&simple_info.fc_pdf[i]);
			}
			fclose(flin);
		}

		//如果采用固定的fc_ball
		if (controller[0].Command_Exist("sits_constant_fcball"))
		{
			simple_info.is_constant_fc_ball = 1;
			sscanf(controller[0].Command("sits_constant_fcball"),"%f",&simple_info.constant_fc_ball);
		}

		//随机游走的fc上下限和步长
		if (controller[0].Command_Exist("sits_fcball_max"))
		{
			sscanf(controller[0].Command("sits_fcball_max"), "%f", &simple_info.fc_max);
		}
		if (controller[0].Command_Exist("sits_fcball_min"))
		{
			sscanf(controller[0].Command("sits_fcball_min"), "%f", &simple_info.fc_min);
		}
		if (controller[0].Command_Exist("sits_fcball_move_length"))
		{
			sscanf(controller[0].Command("sits_fcball_move_length"), "%f", &simple_info.move_length);
		}

		//随机种子
		if (controller[0].Command_Exist("sits_fcball_random_seed"))
		{
			sscanf(controller[0].Command("sits_fcball_random_seed"), "%d", &simple_info.random_seed);
		}
	}
    //classical SITS模式的额外初始化
    else if (info.sits_mode == CLASSICAL_SITS_MODE)
    {
		classical_info.k_numbers = 40;
		if (controller[0].Command_Exist("sits_temperature_numbers"))
		{
			classical_info.k_numbers = atoi(controller[0].Command("sits_temperature_numbers"));
		}
		printf("	SITS Temperature Numbers is %d\n", classical_info.k_numbers);

		float temp0 = 300.0f;
		if (controller[0].Command_Exist("target_temperature"))
		{
			temp0 = atof(controller[0].Command("target_temperature"));
		}
		classical_info.beta0 = 1.0f / CONSTANT_kB / temp0;

		float temph = 2.0f * temp0;
		if (controller[0].Command_Exist("sits_temperature_high"))
		{
			temph = atof(controller[0].Command("sits_temperature_high"));
		}
		printf("	SITS Temperature High Border is %.2f\n", temph);

		float templ = temp0 / 1.2f;
		if (controller[0].Command_Exist("sits_temperature_low"))
		{
			templ = atof(controller[0].Command("sits_temperature_low"));
		}
		printf("	SITS Temperature Low Border is %.2f\n", templ);

		classical_info.energy_multiple = 1.0f;
		if (controller[0].Command_Exist("sits_energy_multiple"))
		{
			classical_info.energy_multiple = atof(controller[0].Command("sits_energy_multiple"));
			printf("	SITS Energy Will Times %f\n", classical_info.energy_multiple);
		}

		classical_info.energy_shift = 0.0f;
		if (controller[0].Command_Exist("sits_energy_shift"))
		{
			classical_info.energy_shift = atof(controller[0].Command("sits_energy_shift"));
			printf("	SITS Energy Will Add %f\n", classical_info.energy_shift);
		}

		classical_info.fb_shift = 0.0f;
		if (controller[0].Command_Exist("sits_fb_shift"))
		{
			classical_info.fb_shift = atof(controller[0].Command("sits_fb_shift"));
			printf("	SITS Fc Ball Will Add %f\n", classical_info.fb_shift);
		}

		classical_info.constant_nk = 0;
		if (controller[0].Command_Exist("sits_constant_nk"))
		{
			classical_info.constant_nk = atoi(controller[0].Command("sits_constant_nk"));
		}

		if (!classical_info.constant_nk)
		{
			classical_info.record_interval = 1;
			if (controller[0].Command_Exist("sits_record_interval"))
			{
				classical_info.record_interval = atoi(controller[0].Command("sits_record_interval"));
			}
			classical_info.update_interval = 100;
			if (controller[0].Command_Exist("sits_update_interval"))
			{
				classical_info.update_interval = atoi(controller[0].Command("sits_update_interval"));
			}
		}

		float *tempf;
		Malloc_Safely((void**)&tempf, sizeof(float)*classical_info.k_numbers);

		Cuda_Malloc_Safely((void**)&classical_info.beta_k, sizeof(float)*classical_info.k_numbers);
		Cuda_Malloc_Safely((void**)&classical_info.NkExpBetakU, sizeof(float)*classical_info.k_numbers);
		Cuda_Malloc_Safely((void**)&classical_info.Nk, sizeof(float)*classical_info.k_numbers);
		Cuda_Malloc_Safely((void**)&classical_info.sum_a, sizeof(float));
		Cuda_Malloc_Safely((void**)&classical_info.sum_b, sizeof(float));
		Cuda_Malloc_Safely((void**)&classical_info.d_fc_ball, sizeof(float) * 2);  //这里分配两个，一个存上一次的一个，免得变化太大体系炸了
		Reset_List << <1, 2 >> >(2, classical_info.d_fc_ball, 1.0);
		Cuda_Malloc_Safely((void**)&classical_info.ene_recorded, sizeof(float));
		Cuda_Malloc_Safely((void**)&classical_info.gf, sizeof(float)* classical_info.k_numbers);
		Cuda_Malloc_Safely((void**)&classical_info.gfsum, sizeof(float));
		Cuda_Malloc_Safely((void**)&classical_info.log_weight, sizeof(float)* classical_info.k_numbers);
		Cuda_Malloc_Safely((void**)&classical_info.log_mk_inverse, sizeof(float)* classical_info.k_numbers);
		Cuda_Malloc_Safely((void**)&classical_info.log_norm_old, sizeof(float)* classical_info.k_numbers);
		Cuda_Malloc_Safely((void**)&classical_info.log_norm, sizeof(float)* classical_info.k_numbers);
		Cuda_Malloc_Safely((void**)&classical_info.log_pk, sizeof(float)* classical_info.k_numbers);
		Cuda_Malloc_Safely((void**)&classical_info.log_nk_inverse, sizeof(float)* classical_info.k_numbers);
		Cuda_Malloc_Safely((void**)&classical_info.log_nk, sizeof(float)* classical_info.k_numbers);
		Malloc_Safely((void**)&classical_info.log_nk_recorded_cpu, sizeof(float)*classical_info.k_numbers);
		Malloc_Safely((void**)&classical_info.log_norm_recorded_cpu, sizeof(float)*classical_info.k_numbers);
        if (!classical_info.constant_nk)
        {
            //nk的轨迹文件
            if (controller[0].Command_Exist("sits_nk_traj_file"))
            {
                printf("	SITS Log Nk Trajectory File: %s\n", controller[0].Command("sits_nk_traj_file"));
                Open_File_Safely(&classical_info.nk_traj_file, controller[0].Command("sits_nk_traj_file"), "wb");
            }
            else
            {
                printf("	SITS Log Nk Trajectory File: sits_nk_traj_file.dat\n");
                Open_File_Safely(&classical_info.nk_traj_file, "sits_nk_traj_file.dat", "wb");
            }
            //norm的轨迹文件
            if (controller[0].Command_Exist("sits_norm_traj_file"))
            {
                printf("	SITS Log Normalization Trajectory File: %s\n", controller[0].Command("sits_norm_traj_file"));
                Open_File_Safely(&classical_info.norm_traj_file, controller[0].Command("sits_norm_traj_file"), "wb");
            }
            else
            {
                printf("	SITS Log Normalization Trajectory File: sits_norm_traj_file.dat\n");
                Open_File_Safely(&classical_info.norm_traj_file, "sits_norm_traj_file.dat", "wb");
            }
        }
		//nk的重开文件
		if (controller[0].Command_Exist("sits_nk_rest_file"))
		{
			printf("	SITS Log Nk Restart File: %s\n", controller[0].Command("sits_nk_rest_file"));
			strcpy(classical_info.nk_rest_file, controller[0].Command("sits_nk_rest_file"));
		}
		else
		{
			printf("	SITS Log Nk Restart File: sits_nk_rest_file.dat\n");
			strcpy(classical_info.nk_rest_file, "sits_nk_rest_file.dat");
		}
		//norm的重开文件
		if (controller[0].Command_Exist("sits_norm_rest_file"))
		{
			printf("	SITS Log Normalization Restart File: %s\n", controller[0].Command("sits_norm_rest_file"));
			strcpy(classical_info.nk_rest_file, controller[0].Command("sits_norm_rest_file"));
		}
		else
		{
			printf("	SITS Log Normalization Restart File: sits_norm_rest_file.dat\n");
			strcpy(classical_info.norm_rest_file, "sits_norm_rest_file.dat");
		}
		//nk的初始化文件及其初始化
		if (controller[0].Command_Exist("sits_nk_init_file"))
		{
			FILE *nk_init_file;
			printf("	SITS Log Nk Initial File: %s\n", controller[0].Command("sits_nk_init_file"));
			Open_File_Safely(&nk_init_file, controller[0].Command("sits_nk_init_file"), "wb");
			for (int i = 0; i < classical_info.k_numbers; i++)
			{
				fscanf(nk_init_file, "%f", &tempf[i]);
			}
			fclose(nk_init_file);
		}
		else
		{
			printf("	SITS Log Nk Initial To Default Value 0.0\n");
			for (int i = 0; i < classical_info.k_numbers; i++)
			{
				tempf[i] = 0.0; 
			}
		}
		cudaMemcpy(classical_info.log_nk, tempf, sizeof(float)* classical_info.k_numbers, cudaMemcpyHostToDevice);
		
		for (int i = 0; i < classical_info.k_numbers; i++)
		{
			tempf[i] = expf(tempf[i]);
		}
		cudaMemcpy(classical_info.Nk, tempf, sizeof(float)* classical_info.k_numbers, cudaMemcpyHostToDevice);

		//norm的初始化文件及其初始化
		if (controller[0].Command_Exist("sits_norm_init_file"))
		{
			FILE *norm_init_file;
			printf("	SITS Log Normalization Initial File: %s\n", controller[0].Command("sits_norm_init_file"));
			Open_File_Safely(&norm_init_file, controller[0].Command("sits_norm_init_file"), "wb");
			for (int i = 0; i < classical_info.k_numbers; i++)
			{
				fscanf(norm_init_file, "%f", &tempf[i]);
			}
			fclose(norm_init_file);
		}
		else
		{
			printf("	SITS Log Normalization Initial To Default Value %.0e\n", -FLT_MAX);
			for (int i = 0; i < classical_info.k_numbers; i++)
			{
				tempf[i] = -FLT_MAX;
			}
		}
		cudaMemcpy(classical_info.log_norm, tempf, sizeof(float)* classical_info.k_numbers, cudaMemcpyHostToDevice);
		cudaMemcpy(classical_info.log_norm_old, tempf, sizeof(float)* classical_info.k_numbers, cudaMemcpyHostToDevice);

		//温度相关信息
		float temp_slope = (temph - templ) / (classical_info.k_numbers - 1);
		for (int i = 0; i < classical_info.k_numbers; i = i + 1)
		{
			tempf[i] = templ + temp_slope * i;
			tempf[i] = 1. / (CONSTANT_kB * tempf[i]);
		}
		cudaMemcpy(classical_info.beta_k, tempf, sizeof(float)* classical_info.k_numbers, cudaMemcpyHostToDevice);
		

		free(tempf);
    }
	printf("End (Selective Integrated Tempering Sampling)\n\n");
}
void SITS::Clear_SITS_Energy()
{
	Reset_List << <ceilf((float)info.atom_numbers / 32), 32 >> >
		(info.protein_atom_numbers, d_protein_water_atom_energy, 0.);//清空用于记录AB相互作用的列表
}
void SITS::Calculate_Total_SITS_Energy(float *d_atom_energy)
{
	Sum_Of_List << <1, 1024 >> >(0, info.protein_atom_numbers, d_atom_energy, d_total_pp_atom_energy);
	Sum_Of_List << <1, 1024 >> >(info.protein_atom_numbers, info.atom_numbers, d_atom_energy, d_total_ww_atom_energy);
	Sum_Of_List << <1, 1024 >> >(0,info.protein_atom_numbers, d_protein_water_atom_energy, d_total_protein_water_atom_energy);
}

void SITS::Print()
{
	if (info.sits_mode == CLASSICAL_SITS_MODE)
	{
		cudaMemcpy(&info.fc_ball, classical_info.d_fc_ball, sizeof(float), cudaMemcpyDeviceToHost);
	}
    cudaMemcpy(&h_total_pp_atom_energy, d_total_pp_atom_energy, sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_total_ww_atom_energy, d_total_ww_atom_energy, sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_total_protein_water_atom_energy, d_total_protein_water_atom_energy, sizeof(float), cudaMemcpyDeviceToHost);
    printf("SITS: ______AA______ ______BB______ ______AB______ fc_ball\n");
    printf("      %14.4f %14.4f %14.4f %7.4f\n", 
        h_total_pp_atom_energy, h_total_ww_atom_energy, h_total_protein_water_atom_energy, info.fc_ball);
    fprintf(sits_ene_record_out, "%f %f %f %f\n",
        h_total_pp_atom_energy, h_total_ww_atom_energy, h_total_protein_water_atom_energy, info.fc_ball);
}

void SITS::Prepare_For_Calculate_Force(int *need_potential, int isPrintStep)
{
	Reset_List << <ceilf((float)3.*info.atom_numbers / 128), 128 >> >(3 * info.atom_numbers, (float*)protein_water_frc, 0.);
	if (info.sits_mode <= 1 && isPrintStep)
		*need_potential += 1;
	if (*need_potential > 0)
		Clear_SITS_Energy();
}

void SITS::SITS_Enhanced_Force(int steps, VECTOR *frc)
{
	if (info.sits_mode == SIMPLE_SITS_MODE)
	{
		//确定fc_ball数值（通过在给定势场中随机游走刷新，以保证fc_ball的制定分布）
		if (!simple_info.is_constant_fc_ball)
		{
			fc_ball_random_walk();
		}
		else
		{
			info.fc_ball = simple_info.constant_fc_ball;
		}
		SITS_For_Enhanced_Force_Protein << <1, 128 >> >
			(info.protein_atom_numbers, frc, protein_water_frc, info.fc_ball, info.pwwp_enhance_factor* info.fc_ball + 1.0 - info.pwwp_enhance_factor);
		SITS_For_Enhanced_Force_Water << <1, 128 >> >
			(info.protein_atom_numbers, info.atom_numbers, frc, protein_water_frc, info.pwwp_enhance_factor* info.fc_ball + 1.0 - info.pwwp_enhance_factor);
	}
	else if (info.sits_mode == CLASSICAL_SITS_MODE)
	{
		SITS_Classical_Update_Info(steps);
		SITS_Classical_Enhanced_Force << <1, 1 >> >(info.atom_numbers, info.protein_atom_numbers, 
			info.pwwp_enhance_factor, frc, protein_water_frc, 
			d_total_pp_atom_energy, d_total_protein_water_atom_energy,
			classical_info.k_numbers, classical_info.NkExpBetakU,
			classical_info.beta_k, classical_info.Nk,
			classical_info.sum_a, classical_info.sum_b, classical_info.d_fc_ball,
			classical_info.beta0, classical_info.energy_multiple, 
			classical_info.energy_shift, classical_info.fb_shift);
	}
	else
	{
		SITS_For_Enhanced_Force_Protein << <1, 128 >> >
			(info.protein_atom_numbers, frc, protein_water_frc, info.fc_ball, info.pwwp_enhance_factor* info.fc_ball + 1.0 - info.pwwp_enhance_factor);
		SITS_For_Enhanced_Force_Water << <1, 128 >> >
			(info.protein_atom_numbers, info.atom_numbers, frc, protein_water_frc, info.pwwp_enhance_factor* info.fc_ball + 1.0 - info.pwwp_enhance_factor);
	}
}
float SITS::FC_BALL_INFORMATION::get_fc_probability(float pos)
{

	//float lin = (float)grid_numbers - (pos - fc_min)*grid_numbers / (fc_max - fc_min);//如果fcball推荐分布的列表是反向放置的，那么就按这条代码计算
	float lin = (float) (pos - fc_min)*grid_numbers / (fc_max - fc_min);

	int a = (int)lin;
	if (a >= 0 && a<grid_numbers)
	{
		return fc_pdf[a];
	}
	else
	{
		return 0.;
	}
}

void SITS::fc_ball_random_walk()
{
	float fc_test = info.fc_ball + simple_info.move_length*((float)rand() / RAND_MAX - 0.5);//尝试更新的fc_ball数值
	float  p2 = simple_info.get_fc_probability(fc_test);//新位置的势能
	float p = p2 / simple_info.current_fc_probability;//新旧位置的势能比,（向势场高处走，最后的分布就正比于势能曲线）
	if (p > ((float)rand() / RAND_MAX))
	{
		info.fc_ball = fc_test;
		simple_info.current_fc_probability = p2;
		if (info.fc_ball >= simple_info.fc_max)
		{
			info.fc_ball = simple_info.fc_max;
		}
		else if (info.fc_ball <= simple_info.fc_min)
		{
			info.fc_ball = simple_info.fc_min;
		}
	}
}

void SITS::CLASSICAL_SITS_INFORMATION::Export_Restart_Information_To_File()
{
	FILE *nk;
	Open_File_Safely(&nk, nk_rest_file, "w");
	cudaMemcpy(log_nk_recorded_cpu, Nk, sizeof(float)*k_numbers, cudaMemcpyDeviceToHost);
	for (int i = 0; i < k_numbers; i++)
	{
		fprintf(nk, "%f\n", log_nk_recorded_cpu[i]);
	}
	fclose(nk);
	
	FILE *norm;
	Open_File_Safely(&norm, norm_rest_file, "w");
	cudaMemcpy(log_norm_recorded_cpu, log_norm, sizeof(float)*k_numbers, cudaMemcpyDeviceToHost);
	for (int i = 0; i < k_numbers; i++)
	{
		fprintf(norm, "%f\n", log_norm_recorded_cpu[i]);
	}
	fclose(norm);

}

void SITS::Clear_SITS()
{
	if (sits_ene_record_out != NULL)
	{
		fclose(sits_ene_record_out);
	}

	if ((info.sits_mode & 0x0000000F) == SIMPLE_SITS_MODE)
	{
		if (simple_info.fc_pdf != NULL)
		{
			free(simple_info.fc_pdf);
		}
	}

	if (d_total_pp_atom_energy != NULL)
	{
		cudaFree(d_total_pp_atom_energy);
	}
	if (d_total_ww_atom_energy != NULL)
	{
		cudaFree(d_total_ww_atom_energy);
	}
	if (d_protein_water_atom_energy != NULL)
	{
		cudaFree(d_protein_water_atom_energy);
	}
	if (d_total_protein_water_atom_energy != NULL)
	{
		cudaFree(d_total_protein_water_atom_energy);
	}

	if (protein_water_frc != NULL)
	{
		cudaFree(protein_water_frc);
	}

	ATOM_GROUP *h_nl = NULL;//临时的近邻表
	Malloc_Safely((void**)&h_nl, sizeof(ATOM_GROUP)*info.atom_numbers);
	if (d_nl_ppww != NULL)
	{
		cudaMemcpy(h_nl, d_nl_ppww, sizeof(ATOM_GROUP)*info.atom_numbers, cudaMemcpyDeviceToHost);
		for (int i = 0; i < info.atom_numbers; i = i + 1)
		{
			cudaFree(h_nl[i].atom_serial);
		}
		cudaFree(d_nl_ppww);
	}
	if (d_nl_pwwp != NULL)
	{
		cudaMemcpy(h_nl, d_nl_ppww, sizeof(ATOM_GROUP)*info.atom_numbers, cudaMemcpyDeviceToHost);
		for (int i = 0; i < info.atom_numbers; i = i + 1)
		{
			cudaFree(h_nl[i].atom_serial);
		}
		cudaFree(d_nl_pwwp);
	}
	free(h_nl);
}
