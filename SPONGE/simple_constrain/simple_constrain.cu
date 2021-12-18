#include "simple_constrain.cuh"
__global__ void Constrain_Force_Cycle
(const int constrain_pair_numbers, const UNSIGNED_INT_VECTOR *uint_crd, const VECTOR scaler,
const CONSTRAIN_PAIR *constrain_pair,const VECTOR *pair_dr,
	VECTOR *test_frc)
{
	int pair_i = blockDim.x*blockIdx.x + threadIdx.x;
	if (pair_i < constrain_pair_numbers)
	{
		CONSTRAIN_PAIR cp = constrain_pair[pair_i];
		float r_1;
		VECTOR dr;
		float frc_abs;
		VECTOR frc_lin;

		dr.x = ((int)(uint_crd[cp.atom_i_serial].uint_x - uint_crd[cp.atom_j_serial].uint_x)) * scaler.x;
		dr.y = ((int)(uint_crd[cp.atom_i_serial].uint_y - uint_crd[cp.atom_j_serial].uint_y)) * scaler.y;
		dr.z = ((int)(uint_crd[cp.atom_i_serial].uint_z - uint_crd[cp.atom_j_serial].uint_z)) * scaler.z;
		r_1=rnorm3df(dr.x, dr.y, dr.z);
		frc_abs = (1. - cp.constant_r*r_1)*cp.constrain_k;

		frc_lin.x = frc_abs*pair_dr[pair_i].x;
		frc_lin.y = frc_abs*pair_dr[pair_i].y;
		frc_lin.z = frc_abs*pair_dr[pair_i].z;


		atomicAdd(&test_frc[cp.atom_j_serial].x, frc_lin.x);
		atomicAdd(&test_frc[cp.atom_j_serial].y, frc_lin.y);
		atomicAdd(&test_frc[cp.atom_j_serial].z, frc_lin.z);

		atomicAdd(&test_frc[cp.atom_i_serial].x, -frc_lin.x);
		atomicAdd(&test_frc[cp.atom_i_serial].y, -frc_lin.y);
		atomicAdd(&test_frc[cp.atom_i_serial].z, -frc_lin.z);
	}
}
__global__ void Refresh_Uint_Crd(const int atom_numbers, const VECTOR *crd, const VECTOR quarter_crd_to_uint_crd_cof, UNSIGNED_INT_VECTOR *uint_crd, const VECTOR *test_frc,
	const float *mass_inverse,const float half_exp_gamma_plus_half)
{
	int atom_i = blockDim.x*blockIdx.x + threadIdx.x;
	if (atom_i < atom_numbers)
	{
		INT_VECTOR tempi;
		VECTOR crd_lin = crd[atom_i];
		VECTOR frc_lin = test_frc[atom_i];
		float mass_lin = mass_inverse[atom_i];

		crd_lin.x = crd_lin.x + half_exp_gamma_plus_half*frc_lin.x*mass_lin;//mass实际为mass的倒数，frc_lin已经乘以dt^2
		crd_lin.y = crd_lin.y + half_exp_gamma_plus_half*frc_lin.y*mass_lin;
		crd_lin.z = crd_lin.z + half_exp_gamma_plus_half*frc_lin.z*mass_lin;

		tempi.int_x = crd_lin.x*quarter_crd_to_uint_crd_cof.x;
		tempi.int_y = crd_lin.y*quarter_crd_to_uint_crd_cof.y;
		tempi.int_z = crd_lin.z*quarter_crd_to_uint_crd_cof.z;

		uint_crd[atom_i].uint_x = tempi.int_x << 2;
		uint_crd[atom_i].uint_y = tempi.int_y << 2;
		uint_crd[atom_i].uint_z = tempi.int_z << 2;
	}
}
__global__ void Last_Crd_To_dr
(const int constarin_pair_numbers, const VECTOR *atom_crd,
const VECTOR quarter_crd_to_uint_crd_cof, const VECTOR uint_dr_to_dr,
const CONSTRAIN_PAIR *constrain_pair,
VECTOR *pair_dr)
{
	int pair_i = blockDim.x*blockIdx.x + threadIdx.x;
	if (pair_i < constarin_pair_numbers)
	{
		INT_VECTOR tempi;
		INT_VECTOR tempj;
		UNSIGNED_INT_VECTOR uint_crd_i;
		UNSIGNED_INT_VECTOR uint_crd_j;
		CONSTRAIN_PAIR cp = constrain_pair[pair_i];
		VECTOR dr;

		tempi.int_x = atom_crd[cp.atom_i_serial].x*quarter_crd_to_uint_crd_cof.x;
		tempi.int_y = atom_crd[cp.atom_i_serial].y*quarter_crd_to_uint_crd_cof.y;
		tempi.int_z = atom_crd[cp.atom_i_serial].z*quarter_crd_to_uint_crd_cof.z;

		tempj.int_x = atom_crd[cp.atom_j_serial].x*quarter_crd_to_uint_crd_cof.x;
		tempj.int_y = atom_crd[cp.atom_j_serial].y*quarter_crd_to_uint_crd_cof.y;
		tempj.int_z = atom_crd[cp.atom_j_serial].z*quarter_crd_to_uint_crd_cof.z;

		uint_crd_i.uint_x = tempi.int_x << 2;
		uint_crd_i.uint_y = tempi.int_y << 2;
		uint_crd_i.uint_z = tempi.int_z << 2;

		uint_crd_j.uint_x = tempj.int_x << 2;
		uint_crd_j.uint_y = tempj.int_y << 2;
		uint_crd_j.uint_z = tempj.int_z << 2;

		dr.x = ((int)(uint_crd_i.uint_x - uint_crd_j.uint_x)) * uint_dr_to_dr.x;
		dr.y = ((int)(uint_crd_i.uint_y - uint_crd_j.uint_y)) * uint_dr_to_dr.y;
		dr.z = ((int)(uint_crd_i.uint_z - uint_crd_j.uint_z)) * uint_dr_to_dr.z;

		pair_dr[pair_i] = dr;
	}
}
__global__ void Refresh_Crd_Vel(const int atom_numbers, const float dt_inverse,const float dt, VECTOR *crd, VECTOR *vel, const VECTOR *test_frc,
	const float *mass_inverse, const float exp_gamma, const float half_exp_gamma_plus_half)
{
	int atom_i = blockDim.x*blockIdx.x + threadIdx.x;
	if (atom_i < atom_numbers)
	{
		VECTOR crd_lin = crd[atom_i];
		VECTOR frc_lin = test_frc[atom_i];
		VECTOR vel_lin = vel[atom_i];
		float mass_lin = mass_inverse[atom_i];

		frc_lin.x = frc_lin.x*mass_lin;
		frc_lin.y = frc_lin.y*mass_lin;
		frc_lin.z = frc_lin.z*mass_lin;//mass实际为mass的倒数，frc_lin已经乘以dt^2

		crd_lin.x = crd_lin.x + half_exp_gamma_plus_half*frc_lin.x;
		crd_lin.y = crd_lin.y + half_exp_gamma_plus_half*frc_lin.y;
		crd_lin.z = crd_lin.z + half_exp_gamma_plus_half*frc_lin.z;


		vel_lin.x = (vel_lin.x + exp_gamma*frc_lin.x*dt_inverse);
		vel_lin.y = (vel_lin.y + exp_gamma*frc_lin.y*dt_inverse);
		vel_lin.z = (vel_lin.z + exp_gamma*frc_lin.z*dt_inverse);

		crd[atom_i] = crd_lin;
		vel[atom_i] = vel_lin;
	}
}

void SIMPLE_CONSTARIN::Initial_Simple_Constrain(CONTROLLER *controller, const int atom_numbers, const float dt, const VECTOR box_length, const float exp_gamma, const int is_Minimization, float *atom_mass, int *system_freedom)
{
	//从传入的参数复制基本信息
	info.atom_numbers = atom_numbers;
	info.dt = dt;
	info.quarter_crd_to_uint_crd_cof = 0.25 * CONSTANT_UINT_MAX_FLOAT / box_length;
	info.uint_dr_to_dr_cof = 1.0f / CONSTANT_UINT_MAX_FLOAT * box_length;
	info.volume = box_length.x * box_length.y * box_length.z;


	//确定使用的朗之万热浴，并给定exp_gamma
	info.exp_gamma = exp_gamma;
	half_exp_gamma_plus_half = 0.5*(1. + info.exp_gamma);
	if (is_Minimization)
	{
		info.exp_gamma = 0.0f;
	}
	
	info.iteration_numbers = 25;
	if (controller[0].Command_Exist(this->module_name, "iteration_numbers"))
	{
		sscanf(controller[0].Command(this->module_name, "iteration_numbers"), "%d", &info.iteration_numbers);
	}
	controller[0].printf("    constrain iteration step is %d\n", info.iteration_numbers);
	if (controller[0].Command_Exist(this->module_name, "step_length"))
	{
		sscanf(controller[0].Command(this->module_name, "step_length"), "%f", &info.step_length);
	}
	controller[0].printf("    constrain step length is %.2f\n", info.step_length);
	
	int extra_numbers = 0;
	FILE *fp = NULL;
	//读文件第一个数确认constrain数量，为分配内存做准备
	if (controller[0].Command_Exist(this->module_name, "in_file"))
	{
		Open_File_Safely(&fp, controller[0].Command(this->module_name, "in_file"), "r");
		fscanf(fp, "%d", &extra_numbers);
	}

	dt_inverse = 1. / info.dt;
	constrain_pair_numbers = bond_constrain_pair_numbers + angle_constrain_pair_numbers + extra_numbers;

	system_freedom[0] -= constrain_pair_numbers;

	controller[0].printf("    constrain pair number is %d\n", constrain_pair_numbers);

	Cuda_Malloc_Safely((void**)&constrain_frc, sizeof(VECTOR)*info.atom_numbers);
	Cuda_Malloc_Safely((void**)&test_uint_crd, sizeof(UNSIGNED_INT_VECTOR)*info.atom_numbers);
	Cuda_Malloc_Safely((void**)&last_pair_dr, sizeof(VECTOR)*constrain_pair_numbers);
	Cuda_Malloc_Safely((void**)&d_pair_virial, sizeof(float)*constrain_pair_numbers);
	Cuda_Malloc_Safely((void**)&d_virial, sizeof(float));

	Malloc_Safely((void**)&h_constrain_pair, sizeof(CONSTRAIN_PAIR)*constrain_pair_numbers);
	Cuda_Malloc_Safely((void**)&constrain_pair, sizeof(CONSTRAIN_PAIR)*constrain_pair_numbers);
	for (int i = 0; i < bond_constrain_pair_numbers; i = i + 1)
	{
		h_constrain_pair[i] = h_bond_pair[i];
		h_constrain_pair[i].constrain_k = info.step_length / half_exp_gamma_plus_half
			* h_constrain_pair[i].constrain_k;
	}
	for (int i = 0; i < angle_constrain_pair_numbers; i = i + 1)
	{
		h_constrain_pair[i + bond_constrain_pair_numbers] = h_angle_pair[i];
		h_constrain_pair[i + bond_constrain_pair_numbers].constrain_k = info.step_length/half_exp_gamma_plus_half
			* h_constrain_pair[i + bond_constrain_pair_numbers].constrain_k;
	}
	//读文件存入
	if (fp != NULL)
	{
		int atom_i, atom_j;
		int count = bond_constrain_pair_numbers + angle_constrain_pair_numbers;
		for (int i = 0; i < extra_numbers; i = i + 1)
		{
			fscanf(fp, "%d %d %f", &atom_i, &atom_j, &h_constrain_pair[count].constant_r);
			h_constrain_pair[count].atom_i_serial = atom_i;
			h_constrain_pair[count].atom_j_serial = atom_j;
			h_constrain_pair[count].constrain_k = info.step_length / half_exp_gamma_plus_half
				* atom_mass[atom_i] * atom_mass[atom_j] / (atom_mass[atom_i] + atom_mass[atom_j]);
			count += 1;
		}
		fclose(fp);
		fp = NULL;
	}

	//上传GPU
	cudaMemcpy(constrain_pair, h_constrain_pair, sizeof(CONSTRAIN_PAIR)*constrain_pair_numbers, cudaMemcpyHostToDevice);


	//清空初始化时使用的临时变量
	if (h_bond_pair != NULL)
	{
		free(h_bond_pair);
		h_bond_pair = NULL;
	}
	if (h_angle_pair != NULL)
	{
		free(h_angle_pair);
		h_angle_pair = NULL;
	}

	controller[0].printf("END INITIALIZING SIMPLE CONSTRAIN\n\n");
	is_initialized = 1;
}

void SIMPLE_CONSTARIN::Remember_Last_Coordinates(VECTOR *crd)
{
	if (is_initialized)
	{
		//获得分子模拟迭代中上一步的距离信息
		Last_Crd_To_dr << <ceilf((float)constrain_pair_numbers / 128), 128 >> >
			(constrain_pair_numbers, crd,
			info.quarter_crd_to_uint_crd_cof, info.uint_dr_to_dr_cof,
			constrain_pair,
			last_pair_dr);
	}
}

__global__ void Constrain_Force_Cycle_With_Virial
(const int constrain_pair_numbers, const UNSIGNED_INT_VECTOR *uint_crd, const VECTOR scaler,
const CONSTRAIN_PAIR *constrain_pair, const VECTOR *pair_dr,
VECTOR *test_frc, float *d_atom_virial)
{
	int pair_i = blockDim.x*blockIdx.x + threadIdx.x;
	if (pair_i < constrain_pair_numbers)
	{
		CONSTRAIN_PAIR cp = constrain_pair[pair_i];
		VECTOR dr0 = pair_dr[pair_i];
		VECTOR dr = Get_Periodic_Displacement(uint_crd[cp.atom_i_serial], uint_crd[cp.atom_j_serial], scaler);
		float r_1 = rnorm3df(dr.x, dr.y, dr.z);
		float frc_abs = (1. - cp.constant_r*r_1)*cp.constrain_k;
		VECTOR frc_lin = frc_abs * dr0;
		d_atom_virial[pair_i] -= frc_lin * dr0;
		//atomicAdd(d_atom_virial + cp.atom_j_serial, -frc_lin * dr0);

		atomicAdd(&test_frc[cp.atom_j_serial].x, frc_lin.x);
		atomicAdd(&test_frc[cp.atom_j_serial].y, frc_lin.y);
		atomicAdd(&test_frc[cp.atom_j_serial].z, frc_lin.z);

		atomicAdd(&test_frc[cp.atom_i_serial].x, -frc_lin.x);
		atomicAdd(&test_frc[cp.atom_i_serial].y, -frc_lin.y);
		atomicAdd(&test_frc[cp.atom_i_serial].z, -frc_lin.z);
	}
}

static __global__ void pressure_fix(float *pressure, float *virial, float factor)
{
	pressure[0] += factor * virial[0];
}

void SIMPLE_CONSTARIN::Constrain
(VECTOR *crd, VECTOR *vel, const float *mass_inverse, int need_pressure, float *d_pressure)
{
	if (is_initialized)
	{
		//清空约束力和维里
		Reset_List << <ceilf((float)3.*info.atom_numbers / 128), 128 >> >
			(3 * info.atom_numbers, (float*)constrain_frc, 0.);
		if (need_pressure > 0)
		{
			Reset_List << <ceilf((float)constrain_pair_numbers / 1024.0f), 1024 >> >(constrain_pair_numbers, d_pair_virial, 0.0f);
			Reset_List << <1, 1 >> >(1, d_virial, 0.0f);
		}
		for (int i = 0; i < info.iteration_numbers; i = i + 1)
		{
		
			Refresh_Uint_Crd << <ceilf((float)info.atom_numbers / 128), 128 >> >
				(info.atom_numbers, crd, info.quarter_crd_to_uint_crd_cof, test_uint_crd, constrain_frc,
				mass_inverse, half_exp_gamma_plus_half);

			if (need_pressure > 0)
			{
				Constrain_Force_Cycle_With_Virial << <ceilf((float)constrain_pair_numbers / 128), 128 >> >
					(constrain_pair_numbers, test_uint_crd, info.uint_dr_to_dr_cof,
					constrain_pair, last_pair_dr,
					constrain_frc, d_pair_virial);
			}
			else
			{
				Constrain_Force_Cycle << <ceilf((float)constrain_pair_numbers / 128), 128 >> >
					(constrain_pair_numbers, test_uint_crd, info.uint_dr_to_dr_cof,
					constrain_pair, last_pair_dr,
					constrain_frc);
			}
		}
		if (need_pressure > 0)
		{
			Sum_Of_List << <1, 1024 >> >(constrain_pair_numbers, d_pair_virial, d_virial);
			pressure_fix << <1, 1 >> >(d_pressure, d_virial, 1 / info.dt / info.dt / 3.0 / info.volume);
		}
	
		Refresh_Crd_Vel << <ceilf((float)info.atom_numbers / 128), 128 >> >
			(info.atom_numbers, dt_inverse, info.dt, crd, vel, constrain_frc,
			mass_inverse, info.exp_gamma, half_exp_gamma_plus_half);
	}
}

void SIMPLE_CONSTARIN::Add_HBond_To_Constrain_Pair
(CONTROLLER *controller, const int bond_numbers, const int *atom_a, const int *atom_b, const float *bond_r,
const float *atom_mass, char *module_name)
{
	controller[0].printf("START INITIALIZING SIMPLE CONSTRAIN:\n");
	if (module_name == NULL)
	{
		strcpy(this->module_name, "simple_constrain");
	}
	else
	{
		strcpy(this->module_name, module_name);
	}
	info.constrain_mass = 3.0f;
	if (controller[0].Command_Exist(this->module_name,"in_file"))
		info.constrain_mass = 0.0f;
	if (controller[0].Command_Exist(this->module_name, "mass"))
		info.constrain_mass = atof(controller[0].Command(this->module_name, "mass"));
	//预先分配一个足够大的CONSTRAIN_PAIR用于临时存储
	Malloc_Safely((void**)&h_bond_pair, sizeof(CONSTRAIN_PAIR)*bond_numbers);
	int s = 0;
	float mass_a, mass_b;
	for (int i = 0; i < bond_numbers; i = i + 1)
	{
		mass_a = atom_mass[atom_a[i]];
		mass_b = atom_mass[atom_b[i]];
		if ((mass_a <info.constrain_mass && mass_a > 0) || (mass_b <info.constrain_mass && mass_b > 0))//含有H原子的bond
		{
			h_bond_pair[s].atom_i_serial = atom_a[i];
			h_bond_pair[s].atom_j_serial = atom_b[i];
			h_bond_pair[s].constant_r = bond_r[i];
			h_bond_pair[s].constrain_k = atom_mass[atom_a[i]] * atom_mass[atom_b[i]] / (atom_mass[atom_a[i]] + atom_mass[atom_b[i]]);
			s = s + 1;
		}
	}
	bond_constrain_pair_numbers = s;

	//假设使用者不会在调用Add_HBond_To_Constrain_Pair与Add_HAngle_To_Constrain_Pair中途释放atom_a等指针指向的空间
	bond_info.bond_numbers = bond_numbers;
	bond_info.atom_a = atom_a;
	bond_info.atom_b = atom_b;
	bond_info.bond_r = bond_r;
}
void SIMPLE_CONSTARIN::Add_HAngle_To_Constrain_Pair
(CONTROLLER *controller, const int angle_numbers, const int *atom_a, const int *atom_b, const int *atom_c,
const float *angle_theta, const float *atom_mass)
{
	//默认认为已经运行了Add_HBond_To_Constrain_Pair

	//预先分配一个足够大的CONSTRAIN_PAIR用于临时存储
	Malloc_Safely((void**)&h_angle_pair, sizeof(CONSTRAIN_PAIR)*angle_numbers*2);
	int s = 0;
	float mass_a, mass_b, mass_c;
	for (int i = 0; i < 0; i = i + 1)
	{
		mass_a = atom_mass[atom_a[i]];
		mass_c = atom_mass[atom_c[i]];
		if ((mass_a <info.constrain_mass && mass_a > 0) || (mass_c <info.constrain_mass && mass_c > 0))//含有H原子的angle,假设氢原子不会在角中心。
		{
			h_angle_pair[s].atom_i_serial = atom_a[i];//固定angle两端的两个点
			h_angle_pair[s].atom_j_serial = atom_c[i];

			float rab=0., rbc=0.;
			for (int j = 0; j < bond_info.bond_numbers; j = j + 1)
			{
				//找到a，b原子的平衡距离
				if ((bond_info.atom_a[j] == atom_a[i] && bond_info.atom_b[j] == atom_b[i])
					|| (bond_info.atom_a[j] == atom_b[i] && bond_info.atom_b[j] == atom_a[i]))
				{
					rab = bond_info.bond_r[j];
				}

				//找到b，c原子的平衡距离
				if ((bond_info.atom_a[j] == atom_c[i] && bond_info.atom_b[j] == atom_b[i])
					|| (bond_info.atom_a[j] == atom_b[i] && bond_info.atom_b[j] == atom_c[i]))
				{
					rbc = bond_info.bond_r[j];
				}
			}
			if (rab == 0. || rbc == 0.)
			{
				controller[0].printf("    Error: Wrong BOND and ANGLE combination!\n");
				getchar();
				continue;
			}

			//运用余弦定理得到平衡的ac长度用于constrain
			h_angle_pair[s].constant_r = sqrtf(rab*rab + rbc*rbc - 2.*rab*rbc*cosf(angle_theta[i]));
			h_angle_pair[s].constrain_k = atom_mass[atom_a[i]] * atom_mass[atom_c[i]] / (atom_mass[atom_a[i]] + atom_mass[atom_c[i]]);

			s = s + 1;

		}
	}
	angle_constrain_pair_numbers = s;
}

void SIMPLE_CONSTARIN::Update_Volume(VECTOR box_length)
{
	if (is_initialized)
	{
		info.quarter_crd_to_uint_crd_cof = 0.25 * CONSTANT_UINT_MAX_FLOAT / box_length;
		info.uint_dr_to_dr_cof = 1.0f / CONSTANT_UINT_MAX_FLOAT * box_length;
		info.volume = box_length.x * box_length.y * box_length.z;
	}
}

void SIMPLE_CONSTARIN::Clear()
{
	if (is_initialized)
	{
		is_initialized = 0;

		cudaFree(constrain_pair);
		constrain_pair = NULL;

		free(h_constrain_pair);
		h_constrain_pair = NULL;

		cudaFree(last_pair_dr);
		last_pair_dr = NULL;

		cudaFree(constrain_frc);
		constrain_frc = NULL;

		cudaFree(test_uint_crd);
		test_uint_crd = NULL;

		cudaFree(d_pair_virial);
		d_pair_virial = NULL;

		cudaFree(d_virial);
		d_virial = NULL;
	}
}
