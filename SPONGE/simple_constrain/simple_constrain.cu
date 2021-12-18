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
		int scanf_ret = sscanf(controller[0].Command(this->module_name, "iteration_numbers"), "%d", &info.iteration_numbers);
	}
	controller[0].printf("    constrain iteration step is %d\n", info.iteration_numbers);
	if (controller[0].Command_Exist(this->module_name, "step_length"))
	{
		int scanf_ret = sscanf(controller[0].Command(this->module_name, "step_length"), "%f", &info.step_length);
	}
	controller[0].printf("    constrain step length is %.2f\n", info.step_length);

	int extra_numbers = 0;
	FILE *fp = NULL;
	//读文件第一个数确认constrain数量，为分配内存做准备
	if (controller[0].Command_Exist(this->module_name, "in_file"))
	{
		Open_File_Safely(&fp, controller[0].Command(this->module_name, "in_file"), "r");
		int scanf_ret = fscanf(fp, "%d", &extra_numbers);
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
		h_constrain_pair[i + bond_constrain_pair_numbers].constrain_k = info.step_length / half_exp_gamma_plus_half
			* h_constrain_pair[i + bond_constrain_pair_numbers].constrain_k;
	}
	//读文件存入
	if (fp != NULL)
	{
		int atom_i, atom_j;
		int count = bond_constrain_pair_numbers + angle_constrain_pair_numbers;
		for (int i = 0; i < extra_numbers; i = i + 1)
		{
			int scanf_ret = fscanf(fp, "%d %d %f", &atom_i, &atom_j, &h_constrain_pair[count].constant_r);
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

	if (controller->Command_Exist(this->module_name, "disable_settle") && atoi(controller->Command(this->module_name, "disable_settle"))) //命令不加settle
	{
		controller->printf("    settle is NOT used for rigid triangle and rigid bond\n");
	}
	else
	{
		controller->printf("    settle is used for rigid triangle and rigid bond\n");
		settle.Initial(controller, this, atom_mass);
	}
	if (is_initialized && !is_controller_printf_initialized)
	{
		is_controller_printf_initialized = 1;
		controller[0].printf("    structure last modify date is %d\n", last_modify_date);
	}
	controller[0].printf("END INITIALIZING SIMPLE CONSTRAIN\n\n");
	is_initialized = 1;
}

void SIMPLE_CONSTARIN::Remember_Last_Coordinates(VECTOR *crd, UNSIGNED_INT_VECTOR *uint_crd, VECTOR scaler)
{
	if (is_initialized)
	{
		//获得分子模拟迭代中上一步的距离信息
		Last_Crd_To_dr << <ceilf((float)constrain_pair_numbers / 128), 128 >> >
			(constrain_pair_numbers, crd,
			info.quarter_crd_to_uint_crd_cof, info.uint_dr_to_dr_cof,
			constrain_pair,
			last_pair_dr);
		settle.Remember_Last_Coordinates(uint_crd, scaler);
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
(VECTOR *crd, VECTOR *vel, const float *mass_inverse, const float *d_mass, VECTOR box_length, int need_pressure, float *d_pressure)
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
		
		settle.Do_SETTLE(d_mass, crd, box_length, vel, need_pressure, d_pressure);
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
	int temp = 0;
	if (controller->Command_Exist(this->module_name, "angle_constrain") && atoi(controller->Command(this->module_name, "angle_constrain")))
	{
		temp = angle_numbers;
	}
	
	//默认认为已经运行了Add_HBond_To_Constrain_Pair

	//预先分配一个足够大的CONSTRAIN_PAIR用于临时存储
	Malloc_Safely((void**)&h_angle_pair, sizeof(CONSTRAIN_PAIR)*angle_numbers*2);
	int s = 0;
	float mass_a, mass_c;
	for (int i = 0; i < temp; i = i + 1)
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

//对几何信息进行转化
//输入：rAB、rAC、rBC：三角形三边长
//输入：mA、mB、mC：ABC三个的质量
//输出：ra rb rc rd re：位置参数，当刚体三角形质心放置于原点时
//A点放置于(0,ra,0)，B点放置于(rc, rb, 0)，C点放置于(rd, re, 0)
__device__ __host__ void Get_Rabcde_From_SSS(float rAB, float rAC, float rBC, float mA, float mB, float mC,
	float &ra, float &rb, float &rc, float &rd, float &re)
{
	float mTotal = mA + mB + mC;
	float Ax = 0;
	float Ay = 0;
	float Bx = -rAB;
	float By = 0;
	float costemp = (rBC * rBC - rAC * rAC - rAB * rAB) / (2 * rAC * rAB);
	float Cx = rAC * costemp;
	float sintemp = sqrtf(1.0f - costemp * costemp);
	float Cy = rAC * sintemp;

	float Ox = (Bx * mB + Cx * mC) / mTotal;
	float Oy = Cy * mC / mTotal;

	Ax -= Ox;
	Ay -= Oy;
	Bx -= Ox;
	By -= Oy;
	Cx -= Ox;
	Cy -= Oy;

	costemp = 1.0f / sqrtf(1.0f + Ax * Ax / Ay / Ay);
	sintemp = costemp * Ax / Ay;

	ra = Ax * sintemp + Ay * costemp;

	rc = Bx * costemp - By * sintemp;
	rb = Bx * sintemp + By * costemp;
	rd = Cx * costemp - Cy * sintemp;
	re = Cx * sintemp + Cy * costemp;

	if (ra < 0)
	{
		ra *= -1;
		rb *= -1;
		re *= -1;
	}
}
//对几何信息进行转化
//输入：rAB、rAC：三角形两边长， angle_BAC：AB和AC的夹角(弧度)
//输入：mA、mB、mC：ABC三个的质量
//输出：ra rb rc rd re：位置参数，当刚体三角形质心放置于原点时
//A点放置于(0,ra,0)，B点放置于(rc, rb, 0)，C点放置于(rd, re, 0)
__device__ __host__ void Get_Rabcde_From_SAS(float rAB, float rAC, float angle_BAC, float mA, float mB, float mC,
	float &ra, float &rb, float &rc, float &rd, float &re)
{
	float mTotal = mA + mB + mC;
	float Ax = 0;
	float Ay = 0;
	float Bx = -rAB;
	float By = 0;

	float costemp = cosf(CONSTANT_Pi - angle_BAC);
	float Cx = rAC * costemp;
	float sintemp = sqrtf(1.0f - costemp * costemp);
	float Cy = rAC * sintemp;

	float Ox = (Bx * mB + Cx * mC) / mTotal;
	float Oy = Cy * mC / mTotal;

	Ax -= Ox;
	Ay -= Oy;
	Bx -= Ox;
	By -= Oy;
	Cx -= Ox;
	Cy -= Oy;

	costemp = 1.0f / sqrtf(1.0f + Ax * Ax / Ay / Ay);
	sintemp = costemp * Ax / Ay;

	ra = Ax * sintemp + Ay * costemp;

	rc = Bx * costemp - By * sintemp;
	rb = Bx * sintemp + By * costemp;
	rd = Cx * costemp - Cy * sintemp;
	re = Cx * sintemp + Cy * costemp;

	if (ra < 0)
	{
		ra *= -1;
		rb *= -1;
		re *= -1;
	}
}

//核心部分
// 部分参考了Shuichi & Peter: SETTLE: An Analytical Version of the SHAKE and RATTLE Algorithm for Rigid Water Models
// A B C 三个点，O 质心
//输入：rB0 上一步的B原子坐标（A为原点）；rC0 上一步的C原子坐标（A为原点）
//rA1 这一步的A原子坐标（质心为原点） rB1 这一步的B原子坐标（质心为原点）；rC1 这一步的C原子坐标（质心为原点）
//ra rb rc rd re：位置参数：当刚体三角形质心放置于原点，A点放置于(0,ra,0)，B点放置于(rc, rb, 0)，C点放置于(rd, re, 0)
// mA、mB、mC：ABC三个的质量 dt:步长
// half_exp_gamma_plus_half, exp_gamma: 同simple_constrain
//输出：rA3 这一步限制后的A原子坐标（质心为原点） rB3 这一步限制后的B原子坐标（质心为 原点）
//rC3 这一步限制后的C原子坐标（质心为原点）vA vB vC 约束后的速度（原位替换）
//virial virial_vector 约束后的维里（原位替换）
__device__ void SETTLE_DO_TRIANGLE(VECTOR rB0, VECTOR rC0,
	VECTOR rA1, VECTOR rB1, VECTOR rC1,
	float ra, float rb, float rc, float rd, float re,
	float mA, float mB, float mC, float dt,
	float half_exp_gamma_plus_half, float exp_gamma,
	VECTOR &rA3, VECTOR &rB3, VECTOR &rC3,
	VECTOR &vA, VECTOR &vB, VECTOR &vC, VECTOR &virial_vector, int triangle_i)
{
	//第0步：构建新坐标系
	//z轴垂直于上一步的BA和BC。 VECTOR ^ VECTOR 是外积
	VECTOR base_vector_z = rB0 ^ rC0;
	//x轴垂直于z轴和这一步的AO
	VECTOR base_vector_x = rA1 ^ base_vector_z;
	//y轴垂直于z轴和x轴
	VECTOR base_vector_y = base_vector_z ^ base_vector_x;
	//归一化
	base_vector_x = rnorm3df(base_vector_x.x, base_vector_x.y, base_vector_x.z) * base_vector_x;
	base_vector_y = rnorm3df(base_vector_y.x, base_vector_y.y, base_vector_y.z) * base_vector_y;
	base_vector_z = rnorm3df(base_vector_z.x, base_vector_z.y, base_vector_z.z) * base_vector_z;

	//第1步：投影至新坐标系
	//     rA0d = {0, 0, 0};
	VECTOR rB0d = { base_vector_x * rB0, base_vector_y * rB0, 0 };
	VECTOR rC0d = { base_vector_x * rC0, base_vector_y * rC0, 0 };
	VECTOR rA1d = { 0, 0, base_vector_z * rA1 };
	VECTOR rB1d = { base_vector_x * rB1, base_vector_y * rB1, base_vector_z * rB1 };
	VECTOR rC1d = { base_vector_x * rC1, base_vector_y * rC1, base_vector_z * rC1 };

	//第2步：绕base_vector_y旋转psi，绕base_vector_x旋转phi得到rX2d
	float sinphi = rA1d.z / ra;
	float cosphi = sqrtf(1.0f - sinphi * sinphi);
	float sinpsi = (rB1d.z - rC1d.z - (rb - re) * sinphi) / ((rd - rc) * cosphi);
	float cospsi = sqrtf(1.0f - sinpsi * sinpsi);

	VECTOR rA2d = { 0.0f, ra * cosphi, rA1d.z };
	VECTOR rB2d = { rc * cospsi, rb * cosphi + rc * sinpsi * sinphi, rB1d.z };
	VECTOR rC2d = { rd * cospsi, re * cosphi + rd * sinpsi * sinphi, rC1d.z };


	//第3步：计算辅助变量 alpha、beta、gamma
	float alpha = rB2d.x * rB0d.x + rC2d.x * rC0d.x + rB2d.y * rB0d.y + rC2d.y * rC0d.y;
	float beta = -rB2d.x * rB0d.y - rC2d.x * rC0d.y + rB2d.y * rB0d.x + rC2d.y * rC0d.x;
	float gamma = rB1d.y * rB0d.x - rB1d.x * rB0d.y + rC1d.y * rC0d.x - rC1d.x * rC0d.y;

	//第4步：绕base_vector_z旋转theta
	float temp = alpha * alpha + beta * beta;
	float sintheta = (alpha * gamma - beta * sqrtf(temp - gamma * gamma)) / temp;
	float costheta = sqrt(1.0f - sintheta * sintheta);
	VECTOR rA3d = { -rA2d.y * sintheta, rA2d.y * costheta, rA2d.z };
	VECTOR rB3d = { rB2d.x * costheta - rB2d.y * sintheta, rB2d.x * sintheta + rB2d.y * costheta, rB2d.z };
	VECTOR rC3d = { rC2d.x * costheta - rC2d.y * sintheta, rC2d.x * sintheta + rC2d.y * costheta, rC2d.z };

	//第5步：投影回去
	rA3 = { rA3d.x * base_vector_x.x + rA3d.y * base_vector_y.x + rA3d.z * base_vector_z.x,
		rA3d.x * base_vector_x.y + rA3d.y * base_vector_y.y + rA3d.z * base_vector_z.y,
		rA3d.x * base_vector_x.z + rA3d.y * base_vector_y.z + rA3d.z * base_vector_z.z };

	rB3 = { rB3d.x * base_vector_x.x + rB3d.y * base_vector_y.x + rB3d.z * base_vector_z.x,
		rB3d.x * base_vector_x.y + rB3d.y * base_vector_y.y + rB3d.z * base_vector_z.y,
		rB3d.x * base_vector_x.z + rB3d.y * base_vector_y.z + rB3d.z * base_vector_z.z };

	rC3 = { rC3d.x * base_vector_x.x + rC3d.y * base_vector_y.x + rC3d.z * base_vector_z.x,
		rC3d.x * base_vector_x.y + rC3d.y * base_vector_y.y + rC3d.z * base_vector_z.y,
		rC3d.x * base_vector_x.z + rC3d.y * base_vector_y.z + rC3d.z * base_vector_z.z };

	//第6步：计算约束造成的速度变化和维里变化
	//节约寄存器，把不用的rX1d拿来当delta vX用
	temp = exp_gamma / dt / half_exp_gamma_plus_half;
	rA1d = temp * (rA3 - rA1);
	rB1d = temp * (rB3 - rB1);
	rC1d = temp * (rC3 - rC1);

	vA = vA + rA1d;
	vB = vB + rB1d;
	vC = vC + rC1d;
	//节约寄存器，把不用的rX0d拿来当FX用
	temp = 1.0f / dt / dt / half_exp_gamma_plus_half;
	//rA0d = temp * mA * (rA3 - rA1);
	rB0d = temp * mB * (rB3 - rB1);
	rC0d = temp * mC * (rC3 - rC1);

	virial_vector.x = rB0d.x * rB0.x + rC0d.x * rC0.x;
	virial_vector.y = rB0d.y * rB0.y + rC0d.y * rC0.y;
	virial_vector.z = rB0d.z * rB0.z + rC0d.z * rC0.z;
}

void SIMPLE_CONSTARIN::SETTLE_INFORMATION::Initial(CONTROLLER* controller, SIMPLE_CONSTARIN* simple_constrain, float *h_mass,
	char* module_name)
{

	if (module_name == NULL)
	{
		strcpy(this->module_name, "settle");
	}
	else
	{
		strcpy(this->module_name, module_name);
	}

	dt = simple_constrain->info.dt;
	half_exp_gamma_plus_half = simple_constrain->half_exp_gamma_plus_half;
	exp_gamma = simple_constrain->info.exp_gamma;

	//遍历搜出simple_constrain里的三角形
	int *linker_numbers = NULL;
	int *linker_atoms = NULL;
	float *link_r = NULL;
	Malloc_Safely((void**)&linker_numbers, sizeof(int)*simple_constrain->info.atom_numbers);
	Malloc_Safely((void**)&linker_atoms, 2 * sizeof(int)*simple_constrain->info.atom_numbers);
	Malloc_Safely((void**)&link_r, 3 * sizeof(float)*simple_constrain->info.atom_numbers);
	for (int i = 0; i < simple_constrain->info.atom_numbers; i++)
	{
		linker_numbers[i] = 0;
	}
	int atom_i, atom_j;
	CONSTRAIN_PAIR pair;
	for (int i = 0; i < simple_constrain->constrain_pair_numbers; i++)
	{
		pair = simple_constrain->h_constrain_pair[i];
		atom_i = pair.atom_i_serial;
		atom_j = pair.atom_j_serial;

		if (linker_numbers[atom_i] < 2 && linker_numbers[atom_j] < 2)
		{
			linker_atoms[2 * atom_i + linker_numbers[atom_i]] = atom_j;
			linker_atoms[2 * atom_j + linker_numbers[atom_j]] = atom_i;
			link_r[3 * atom_i + linker_numbers[atom_i]] = pair.constant_r;
			link_r[3 * atom_j + linker_numbers[atom_j]] = pair.constant_r;
			linker_numbers[atom_i]++;
			linker_numbers[atom_j]++;
		}
		else
		{
			linker_numbers[atom_i] = 3;
			linker_numbers[atom_j] = 3;
		}
	}

	triangle_numbers = 0;
	pair_numbers = 0;
	for (int i = 0; i < simple_constrain->info.atom_numbers; i++)
	{
		if (linker_numbers[i] == 2)
		{
			atom_i = linker_atoms[2 * i];
			atom_j = linker_atoms[2 * i + 1];
			if (linker_numbers[atom_i] == 2 && linker_numbers[atom_j] == 2 &&
				((linker_atoms[2 * atom_i] == i && linker_atoms[2 * atom_i + 1] == atom_j)
				|| (linker_atoms[2 * atom_i + 1] == i && linker_atoms[2 * atom_i] == atom_j)))
			{
				triangle_numbers++;
				linker_numbers[atom_i] = -2;
				linker_numbers[atom_j] = -2;
				if (linker_atoms[2 * atom_i + 1] == atom_j)
				{
					link_r[3 * i + 2] = link_r[3 * atom_i + 1];
				}
				else
				{
					link_r[3 * i + 2] = link_r[3 * atom_i];
				}
			}
			else
			{
				linker_numbers[i] = 3;
				linker_numbers[atom_i] = 3;
				linker_numbers[atom_j] = 3;
			}
		}
		else if (linker_numbers[i] == 1)
		{
			atom_i = linker_atoms[2 * i];
			if (linker_numbers[atom_i] == 1)
			{
				pair_numbers++;
				linker_numbers[atom_i] = -1;
			}
			else
			{
				linker_numbers[i] = 3;
				linker_numbers[atom_i] = 3;
			}
		}
	}

	controller->printf("    rigid triangle numbers is %d\n", triangle_numbers);
	controller->printf("    rigid pair numbers is %d\n", pair_numbers);

	Malloc_Safely((void**)&h_triangles, sizeof(CONSTRAIN_TRIANGLE)* triangle_numbers);
	Cuda_Malloc_Safely((void**)&d_triangles, sizeof(CONSTRAIN_TRIANGLE)* triangle_numbers);
	Malloc_Safely((void**)&h_pairs, sizeof(CONSTRAIN_PAIR)* pair_numbers);
	Cuda_Malloc_Safely((void**)&d_pairs, sizeof(CONSTRAIN_PAIR)* pair_numbers);

	Cuda_Malloc_Safely((void**)&last_triangle_BA, sizeof(VECTOR)* triangle_numbers);
	Cuda_Malloc_Safely((void**)&last_triangle_CA, sizeof(VECTOR)* triangle_numbers);
	Cuda_Malloc_Safely((void**)&last_pair_AB, sizeof(VECTOR)* pair_numbers);
	Cuda_Malloc_Safely((void**)&virial, sizeof(float));
	Cuda_Malloc_Safely((void**)&virial_vector, sizeof(VECTOR)* (triangle_numbers + pair_numbers));
	int triangle_i = 0;
	int pair_i = 0;
	for (int i = 0; i < simple_constrain->info.atom_numbers; i++)
	{
		if (linker_numbers[i] == 2)
		{
			linker_numbers[i] = -2;
			atom_i = linker_atoms[2 * i];
			atom_j = linker_atoms[2 * i + 1];
			h_triangles[triangle_i].atom_A = i;
			h_triangles[triangle_i].atom_B = atom_i;
			h_triangles[triangle_i].atom_C = atom_j;
			Get_Rabcde_From_SSS(link_r[3 * i], link_r[3 * i + 1], link_r[3 * i + 2],
				h_mass[i], h_mass[atom_i], h_mass[atom_j], h_triangles[triangle_i].ra, h_triangles[triangle_i].rb,
				h_triangles[triangle_i].rc, h_triangles[triangle_i].rd, h_triangles[triangle_i].re);
			//printf("%d %d %d %f %f %f\n", i, linker_atoms[2 * i], linker_atoms[2 * i + 1],
			//	link_r[3 * i], link_r[3 * i + 1], link_r[3 * i + 2]);
			triangle_i++;
		}
		if (linker_numbers[i] == 1)
		{
			atom_j = linker_atoms[2 * i];
			linker_numbers[i] = -1;
			h_pairs[pair_i].atom_i_serial = i;
			h_pairs[pair_i].atom_j_serial = atom_j;
			h_pairs[pair_i].constant_r = link_r[3 * i];
			h_pairs[pair_i].constrain_k = 1.0f / (h_mass[i] + h_mass[atom_j]);
			pair_i++;
		}
	}


	cudaMemcpy(d_triangles, h_triangles, sizeof(CONSTRAIN_TRIANGLE)* triangle_numbers, cudaMemcpyHostToDevice);
	cudaMemcpy(d_pairs, h_pairs, sizeof(CONSTRAIN_PAIR)* pair_numbers, cudaMemcpyHostToDevice);

	//原来的重塑
	int new_constrain_pair_numbers = simple_constrain->constrain_pair_numbers - 3 * triangle_numbers - pair_numbers;
	int new_pair_i = 0;

	CONSTRAIN_PAIR *new_h_constrain_pair = NULL;
	Malloc_Safely((void**)&new_h_constrain_pair, sizeof(CONSTRAIN_PAIR)* new_constrain_pair_numbers);

	for (int i = 0; i < simple_constrain->constrain_pair_numbers; i++)
	{
		pair = simple_constrain->h_constrain_pair[i];
		atom_i = pair.atom_i_serial;;
		if (linker_numbers[atom_i] > 0)
		{
			new_h_constrain_pair[new_pair_i] = pair;
			new_pair_i++;
		}
	}
	simple_constrain->constrain_pair_numbers = new_constrain_pair_numbers;
	free(simple_constrain->h_constrain_pair);
	cudaFree(simple_constrain->constrain_pair);
	simple_constrain->h_constrain_pair = new_h_constrain_pair;
	Cuda_Malloc_Safely((void**)&simple_constrain->constrain_pair, sizeof(CONSTRAIN_PAIR)* new_constrain_pair_numbers);
	cudaMemcpy(simple_constrain->constrain_pair, simple_constrain->h_constrain_pair, sizeof(CONSTRAIN_PAIR)* new_constrain_pair_numbers, cudaMemcpyHostToDevice);

	controller->printf("    remaining simple constrain pair numbers is %d\n", new_pair_i);
	for (int i = 0; i < simple_constrain->constrain_pair_numbers; i++)
	{
		pair = simple_constrain->h_constrain_pair[i];
		atom_i = pair.atom_i_serial;
	}
	free(linker_numbers);
	free(linker_atoms);
	free(link_r);
	is_initialized = 1;

}

__global__ void remember_triangle_BA_CA(int triangle_numbers, CONSTRAIN_TRIANGLE* triangles, UNSIGNED_INT_VECTOR *uint_crd, VECTOR scaler,
	VECTOR* last_triangle_BA, VECTOR* last_triangle_CA)
{
	CONSTRAIN_TRIANGLE triangle;
	for (int triangle_i = blockIdx.x * blockDim.x + threadIdx.x; triangle_i < triangle_numbers; triangle_i += blockDim.x * gridDim.x)
	{
		triangle = triangles[triangle_i];
		last_triangle_BA[triangle_i] = Get_Periodic_Displacement(uint_crd[triangle.atom_B], uint_crd[triangle.atom_A], scaler);
		last_triangle_CA[triangle_i] = Get_Periodic_Displacement(uint_crd[triangle.atom_C], uint_crd[triangle.atom_A], scaler);
	}
}

__global__ void remember_pair_AB(int pair_numbers, CONSTRAIN_PAIR* pairs, UNSIGNED_INT_VECTOR *uint_crd, VECTOR scaler,
	VECTOR* last_pair_AB)
{
	CONSTRAIN_PAIR pair;
	for (int pair_i = blockIdx.x * blockDim.x + threadIdx.x; pair_i < pair_numbers; pair_i += blockDim.x * gridDim.x)
	{
		pair = pairs[pair_i];
		last_pair_AB[pair_i] = Get_Periodic_Displacement(uint_crd[pair.atom_j_serial], uint_crd[pair.atom_i_serial], scaler);
	}
}

__global__ void settle_triangle(int triangle_numbers, CONSTRAIN_TRIANGLE* triangles, const float *d_mass,
	VECTOR *crd, VECTOR box_length, VECTOR* last_triangle_BA, VECTOR* last_triangle_CA,
	float dt, float exp_gamma, float half_exp_gamma_plus_half, VECTOR *vel, VECTOR *virial_vector)
{
	CONSTRAIN_TRIANGLE triangle;
	VECTOR rO;
	VECTOR rA, rB, rC;
	float mA, mB, mC;
	for (int triangle_i = blockIdx.x * blockDim.x + threadIdx.x; triangle_i < triangle_numbers; triangle_i += blockDim.x * gridDim.x)
	{
		triangle = triangles[triangle_i];
		rA = crd[triangle.atom_A];
		rB = Get_Periodic_Displacement(crd[triangle.atom_B], rA, box_length);
		rC = Get_Periodic_Displacement(crd[triangle.atom_C], rA, box_length);
		mA = d_mass[triangle.atom_A];
		mB = d_mass[triangle.atom_B];
		mC = d_mass[triangle.atom_C];

		rO = 1.0f / (mA + mB + mC) * (mB * rB + mC * rC) + rA;
		rA = rA - rO;
		rB = rB + rA;
		rC = rC + rA;

		SETTLE_DO_TRIANGLE(last_triangle_BA[triangle_i], last_triangle_CA[triangle_i], rA, rB, rC,
			triangle.ra, triangle.rb, triangle.rc, triangle.rd, triangle.re, mA, mB, mC,
			dt, half_exp_gamma_plus_half, exp_gamma, rA, rB, rC,
			vel[triangle.atom_A], vel[triangle.atom_B], vel[triangle.atom_C], virial_vector[triangle_i], triangle_i);

		crd[triangle.atom_A] = rA + rO;
		crd[triangle.atom_B] = rB + rO;
		crd[triangle.atom_C] = rC + rO;
	}
}

__global__ void settle_pair(int pair_numbers, CONSTRAIN_PAIR* pairs, const float *d_mass,
	VECTOR *crd, VECTOR box_length, VECTOR* last_pair_AB,
	float dt, float exp_gamma, float half_exp_gamma_plus_half,
	VECTOR *vel, VECTOR *virial_vector)
{
	CONSTRAIN_PAIR pair;
	VECTOR r1, r2, kr2;
	float mA, mB, r0r0, r1r1, r1r2, r2r2, k;
	for (int pair_i = blockIdx.x * blockDim.x + threadIdx.x; pair_i < pair_numbers; pair_i += blockDim.x * gridDim.x)
	{
		pair = pairs[pair_i];

		r1 = Get_Periodic_Displacement(crd[pair.atom_j_serial], crd[pair.atom_i_serial], box_length);
		r2 = last_pair_AB[pair_i];
		mA = d_mass[pair.atom_i_serial];
		mB = d_mass[pair.atom_j_serial];

		r0r0 = pair.constant_r * pair.constant_r;
		r1r1 = r1 * r1;
		r1r2 = r1 * r2;
		r2r2 = r2 * r2;

		k = (sqrt(r1r2 * r1r2 - r1r1 * r2r2 + r2r2 * r0r0) - r1r2) / r2r2;
		kr2 = k * r2;

		r1 = -mB * pair.constrain_k * kr2;
		kr2 = mA * pair.constrain_k * kr2;

		crd[pair.atom_i_serial] = crd[pair.atom_i_serial] + r1;
		crd[pair.atom_j_serial] = crd[pair.atom_j_serial] + kr2;

		k = exp_gamma / dt / half_exp_gamma_plus_half;
		vel[pair.atom_i_serial] = vel[pair.atom_i_serial] + k * r1;
		vel[pair.atom_j_serial] = vel[pair.atom_j_serial] + k * kr2;

		r1 = k * mB / dt / dt / half_exp_gamma_plus_half * kr2;
		virial_vector[pair_i].x = r1.x * r2.x;
		virial_vector[pair_i].y = r1.y * r2.y;
		virial_vector[pair_i].z = r1.z * r2.z;
	}
}

void SIMPLE_CONSTARIN::SETTLE_INFORMATION::Remember_Last_Coordinates(UNSIGNED_INT_VECTOR *uint_crd, VECTOR scaler)
{
	if (!is_initialized)
	{
		return;
	}

	remember_pair_AB << <32, 64 >> >(pair_numbers, d_pairs, uint_crd, scaler,
		last_pair_AB);

	remember_triangle_BA_CA << <32, 64 >> >(triangle_numbers, d_triangles, uint_crd, scaler,
		last_triangle_BA, last_triangle_CA);
}

__global__ void Sum_Of_Virial_Vector_To_Pressure(int N, VECTOR* virial_vector, float *pressure, float factor)
{
	float virial = 0;
	for (int i = blockDim.x * blockIdx.x + threadIdx.x; i < N; i = i + blockDim.x * gridDim.x)
	{
		virial = virial + virial_vector[i].x + virial_vector[i].y + virial_vector[i].z;
	}
	atomicAdd(pressure, virial * factor);
}

void SIMPLE_CONSTARIN::SETTLE_INFORMATION::Do_SETTLE(const float *d_mass, VECTOR *crd, VECTOR box_length, VECTOR *vel,
	int need_pressure, float* d_pressure)
{
	if (!is_initialized)
	{
		return;
	}
	settle_pair << <64, 64 >> > (pair_numbers, d_pairs, d_mass, crd, box_length, last_pair_AB,
		dt, exp_gamma, half_exp_gamma_plus_half, vel, virial_vector + triangle_numbers);
	settle_triangle << <64, 64 >> > (triangle_numbers, d_triangles, d_mass, crd, box_length, last_triangle_BA, last_triangle_CA,
		dt, exp_gamma, half_exp_gamma_plus_half, vel, virial_vector);

	if (need_pressure)
	{
		Sum_Of_Virial_Vector_To_Pressure << <64, 1024 >> > (triangle_numbers + pair_numbers, virial_vector, d_pressure, 0.33333f / box_length.x / box_length.y / box_length.z);
	}
}
