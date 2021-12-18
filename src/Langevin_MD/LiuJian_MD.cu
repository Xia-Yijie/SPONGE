#include "LiuJian_MD.cuh"

static __global__ void MD_Iteration_Leap_Frog_With_LiuJian
(const int atom_numbers, const float half_dt, const float dt,const float exp_gamma, 
const float *inverse_mass, const float *sqrt_mass_inverse, 
VECTOR *vel, VECTOR *crd, VECTOR *frc, VECTOR *acc, VECTOR *random_frc)
{
	int i = blockDim.x*blockIdx.x + threadIdx.x;
	if (i < atom_numbers)
	{
		acc[i].x = inverse_mass[i] * frc[i].x;
		acc[i].y = inverse_mass[i] * frc[i].y;
		acc[i].z = inverse_mass[i] * frc[i].z;

		vel[i].x = vel[i].x + dt*acc[i].x;
		vel[i].y = vel[i].y + dt*acc[i].y;
		vel[i].z = vel[i].z + dt*acc[i].z;

		
		crd[i].x = crd[i].x + half_dt*vel[i].x;
		crd[i].y = crd[i].y + half_dt*vel[i].y;
		crd[i].z = crd[i].z + half_dt*vel[i].z;


		vel[i].x = exp_gamma*vel[i].x + sqrt_mass_inverse[i] * random_frc[i].x;
		vel[i].y = exp_gamma*vel[i].y + sqrt_mass_inverse[i] * random_frc[i].y;
		vel[i].z = exp_gamma*vel[i].z + sqrt_mass_inverse[i] * random_frc[i].z;

		crd[i].x = crd[i].x + half_dt*vel[i].x;
		crd[i].y = crd[i].y + half_dt*vel[i].y;
		crd[i].z = crd[i].z + half_dt*vel[i].z;

		frc[i].x = 0.;
		frc[i].y = 0.;
		frc[i].z = 0.;

	}
}
static __global__ void MD_Iteration_Leap_Frog_With_LiuJian_With_Max_Velocity
(const int atom_numbers, const float half_dt, const float dt, const float exp_gamma,
const float *inverse_mass, const float *sqrt_mass_inverse,
VECTOR *vel, VECTOR *crd, VECTOR *frc, VECTOR *acc, VECTOR *random_frc,const float max_vel)
{
	int i = blockDim.x*blockIdx.x + threadIdx.x;
	float abs_vel;
	if (i < atom_numbers)
	{
		acc[i].x = inverse_mass[i] * frc[i].x;
		acc[i].y = inverse_mass[i] * frc[i].y;
		acc[i].z = inverse_mass[i] * frc[i].z;

		vel[i].x = vel[i].x + dt*acc[i].x;
		vel[i].y = vel[i].y + dt*acc[i].y;
		vel[i].z = vel[i].z + dt*acc[i].z;

		abs_vel = norm3df(vel[i].x, vel[i].y, vel[i].z);
		if (abs_vel<max_vel)
		{

		}
		else
		{
			abs_vel = max_vel / abs_vel;
			vel[i].x = abs_vel* vel[i].x;
			vel[i].y = abs_vel* vel[i].y;
			vel[i].z = abs_vel* vel[i].z;
		}

		crd[i].x = crd[i].x + half_dt*vel[i].x;
		crd[i].y = crd[i].y + half_dt*vel[i].y;
		crd[i].z = crd[i].z + half_dt*vel[i].z;


		vel[i].x = exp_gamma*vel[i].x + sqrt_mass_inverse[i] * random_frc[i].x;
		vel[i].y = exp_gamma*vel[i].y + sqrt_mass_inverse[i] * random_frc[i].y;
		vel[i].z = exp_gamma*vel[i].z + sqrt_mass_inverse[i] * random_frc[i].z;

		crd[i].x = crd[i].x + half_dt*vel[i].x;
		crd[i].y = crd[i].y + half_dt*vel[i].y;
		crd[i].z = crd[i].z + half_dt*vel[i].z;

		frc[i].x = 0.;
		frc[i].y = 0.;
		frc[i].z = 0.;

	}
}
void LIUJIAN_MD_INFORMATION::Initial(CONTROLLER *controller, float target_temperature, const char *module_name)
{
	controller[0].printf("START INITIALIZING LANGEVIN_LIU DYNAMICS:\n");
	if (module_name == NULL)
	{
		strcpy(this->module_name, "langevin_liu");
	}
	else
	{
		strcpy(this->module_name, module_name);
	}

	float *h_mass;
	if (!module_name && controller[0].Command_Exist("mass_in_file"))
	{
		FILE *fp;
		Open_File_Safely(&fp, controller[0].Command("mass_in_file"), "r");
		fscanf(fp, "%d", &atom_numbers);
		Malloc_Safely((void**)&h_mass, sizeof(float)*atom_numbers);
		controller[0].printf("    atom_numbers is %d\n", atom_numbers);
		for (int i = 0; i < atom_numbers; i++)
		{
			fscanf(fp, "%f", &h_mass[i]);
		}
		fclose(fp);
	}
	else if (module_name && controller[0].Command_Exist(module_name, "mass_in_file"))
	{
		FILE *fp;
		Open_File_Safely(&fp, controller[0].Command(module_name, "mass_in_file"), "r");
		fscanf(fp, "%d", &atom_numbers);
		Malloc_Safely((void**)&h_mass, sizeof(float)*atom_numbers);
		controller[0].printf("    atom_numbers is %d\n", atom_numbers);
		for (int i = 0; i < atom_numbers; i++)
		{
			fscanf(fp, "%f", &h_mass[i]);
		}
		fclose(fp);
	}
	else if (controller[0].Command_Exist("amber_parm7"))
	{
		controller[0].printf("    Reading mass information from AMBER file:\n");
		Initial_From_AMBER(&h_mass, controller[0].Command("amber_parm7"), controller[0]);
		controller[0].printf("    End reading mass information from AMBER file\n");
	}
	else
	{
		controller[0].printf("    MD basic information needed. Specify the mass in file.\n");
		getchar();
		exit(1);
	}

	this->target_temperature = target_temperature;

	gamma_ln = 1.0f;
	if (controller[0].Command_Exist(this->module_name, "gamma"))
	{
		gamma_ln = atof(controller[0].Command(this->module_name, "gamma"));
	}

	int random_seed = rand();
	if (controller[0].Command_Exist(this->module_name, "seed"))
	{
		random_seed = atoi(controller[0].Command(this->module_name, "seed"));
	}

	controller[0].printf("    target temperature is %.2f K\n", target_temperature);
	controller[0].printf("    friction coefficient is %.2f ps^-1\n", gamma_ln);
	controller[0].printf("    random seed is %d\n", random_seed);

	dt = 0.001;
	if (controller[0].Command_Exist("dt"))
		dt = atof(controller[0].Command("dt"));
	dt *= CONSTANT_TIME_CONVERTION;
	half_dt = 0.5 * dt;

	float4_numbers = ceil((double)3.* atom_numbers / 4.);
	Cuda_Malloc_Safely((void**)&random_force, sizeof(float4)* float4_numbers);
	Cuda_Malloc_Safely((void**)&rand_state, sizeof(curandStatePhilox4_32_10_t)* float4_numbers);

	Setup_Rand_Normal_Kernel << <(unsigned int)ceilf((float)float4_numbers / threads_per_block), threads_per_block >> >
		(float4_numbers, rand_state, random_seed);

	gamma_ln = gamma_ln / CONSTANT_TIME_CONVERTION;//单位换算

	exp_gamma = expf(-gamma_ln * dt);

	float sart_gamma = sqrtf((1. - exp_gamma * exp_gamma) * target_temperature * CONSTANT_kB);
	Cuda_Malloc_Safely((void**)&d_sqrt_mass, sizeof(float)* atom_numbers);
	Cuda_Malloc_Safely((void**)&d_mass_inverse, sizeof(float)* atom_numbers);
	Malloc_Safely((void**)&h_sqrt_mass, sizeof(float)* atom_numbers);
	for (int i = 0; i < atom_numbers; i = i + 1)
	{
		if (h_mass[i] == 0)
			h_sqrt_mass[i] = 0;
		else
			h_sqrt_mass[i] = sart_gamma * sqrtf(1. / h_mass[i]);
	}
	cudaMemcpy(d_sqrt_mass, h_sqrt_mass, sizeof(float)* atom_numbers, cudaMemcpyHostToDevice);

	//确定是否加上速度上限
	max_velocity = 0;
	if (controller[0].Command_Exist(this->module_name, "velocity_max"))
	{
		sscanf(controller[0].Command(this->module_name, "velocity_max"), "%f", &max_velocity);
		controller[0].printf("    max velocity is %.2f\n", max_velocity);
	}
	//记录质量的倒数
	for (int i = 0; i < atom_numbers; i = i + 1)
	{
		if (h_mass[i] == 0)
			h_mass[i] = 0;
		else
			h_mass[i] = 1.0f / h_mass[i];
	}
	cudaMemcpy(d_mass_inverse, h_mass, sizeof(float)*atom_numbers, cudaMemcpyHostToDevice);
	free(h_mass);
	is_initialized = 1;
	if (is_initialized && !is_controller_printf_initialized)
	{
		is_controller_printf_initialized = 1;
		controller[0].printf("    structure last modify date is %d\n", last_modify_date);
	}
	controller[0].printf("END INITIALIZING LANGEVIN_LIU DYNAMICS\n\n");
}

void LIUJIAN_MD_INFORMATION::Initial_From_AMBER(float **h_mass, const char *file_name, CONTROLLER controller)
{
	FILE *parm;
	Open_File_Safely(&parm, file_name, "r");
	char temps[CHAR_LENGTH_MAX];
	char temp_first_str[CHAR_LENGTH_MAX];
	char temp_second_str[CHAR_LENGTH_MAX];
	while (true)
	{
		if (fgets(temps, CHAR_LENGTH_MAX, parm) == NULL)
		{
			break;
		}
		if (sscanf(temps, "%s %s", temp_first_str, temp_second_str) != 2)
		{
			continue;
		}
		//read in atomnumber atomljtypenumber
		if (strcmp(temp_first_str, "%FLAG") == 0
			&& strcmp(temp_second_str, "POINTERS") == 0)
		{
			fgets(temps, CHAR_LENGTH_MAX, parm);

			fscanf(parm, "%d\n", &atom_numbers);
			controller.printf("        atom_numbers is %d\n", atom_numbers);
			Malloc_Safely((void**)h_mass, sizeof(float)*atom_numbers);
		}//FLAG POINTERS

		//atom mass read in
		if (strcmp(temp_first_str, "%FLAG") == 0
			&& strcmp(temp_second_str, "MASS") == 0)
		{
			fgets(temps, CHAR_LENGTH_MAX, parm);
			double lin;
			for (int i = 0; i< atom_numbers; i = i + 1)
			{
				fscanf(parm, "%lf\n", &lin);
				h_mass[0][i] = (float)lin;
			}
		}
	}//while cycle
	fclose(parm);
}

void LIUJIAN_MD_INFORMATION::MD_Iteration_Leap_Frog(VECTOR *frc, VECTOR *vel, VECTOR *acc, VECTOR *crd)
{
	if (is_initialized)
	{
		Rand_Normal << <ceilf((float)float4_numbers / 32.), 32 >> >
			(float4_numbers, rand_state, (float4 *)random_force);

		if (max_velocity <= 0)
		{
			MD_Iteration_Leap_Frog_With_LiuJian
				<< <ceilf((float)atom_numbers / 32), 32 >> >
				(atom_numbers, half_dt, dt, exp_gamma, d_mass_inverse,
				d_sqrt_mass, vel, crd, frc, acc, random_force);
		}
		else
		{
			MD_Iteration_Leap_Frog_With_LiuJian_With_Max_Velocity
				<< <ceilf((float)atom_numbers / 32), 32 >> >
				(atom_numbers, half_dt, dt, exp_gamma, d_mass_inverse,
				d_sqrt_mass, vel, crd, frc, acc, random_force,max_velocity);
			//cudaDeviceSynchronize();
		}
	}
}

void LIUJIAN_MD_INFORMATION::Clear()
{
	if (is_initialized)
	{
		is_initialized = 0;
		cudaFree(rand_state);
		cudaFree(random_force);
		free(h_sqrt_mass);
		cudaFree(d_sqrt_mass);
		cudaFree(d_mass_inverse);

		rand_state = NULL;
		random_force = NULL;
		h_sqrt_mass = NULL;
		d_sqrt_mass = NULL;
		d_mass_inverse = NULL;
	}
}