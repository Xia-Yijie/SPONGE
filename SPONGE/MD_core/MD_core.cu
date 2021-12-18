#include "MD_core.cuh"

#define BOX_TRAJ_COMMAND "box"
#define BOX_TRAJ_DEFAULT_FILENAME "mdbox.txt"
#define TRAJ_COMMAND "crd"
#define TRAJ_DEFAULT_FILENAME "mdcrd.dat"
#define RESTART_COMMAND "rst"
#define RESTART_DEFAULT_FILENAME "restart"

static __global__ void MD_Iteration_Leap_Frog
(const int atom_numbers, VECTOR *vel, VECTOR *crd, VECTOR *frc, VECTOR *acc, const float *inverse_mass, const float dt)
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

		crd[i].x = crd[i].x + dt*vel[i].x;
		crd[i].y = crd[i].y + dt*vel[i].y;
		crd[i].z = crd[i].z + dt*vel[i].z;

		frc[i].x = 0.;
		frc[i].y = 0.;
		frc[i].z = 0.;
	}
}

static __global__ void MD_Iteration_Leap_Frog_With_Max_Velocity
(const int atom_numbers, VECTOR *vel, VECTOR *crd, VECTOR *frc, VECTOR *acc, const float *inverse_mass, const float dt, const float max_velocity)
{
	int i = blockDim.x*blockIdx.x + threadIdx.x;
	if (i < atom_numbers)
	{
		VECTOR acc_i = inverse_mass[i] * frc[i];
		VECTOR vel_i = vel[i] + dt * acc_i;
		vel_i = Make_Vector_Not_Exceed_Value(vel_i, max_velocity);
		vel[i] = vel_i;
		crd[i] = crd[i] + dt * vel_i;
		frc[i] = { 0.0f, 0.0f, 0.0f };
	}
}

static __global__ void MD_Iteration_Gradient_Descent
(const int atom_numbers, VECTOR *crd, VECTOR *frc, const float learning_rate)
{
	int i = blockDim.x*blockIdx.x + threadIdx.x;
	if (i < atom_numbers)
	{
		crd[i].x = crd[i].x + learning_rate * frc[i].x;
		crd[i].y = crd[i].y + learning_rate * frc[i].y;
		crd[i].z = crd[i].z + learning_rate * frc[i].z;

		frc[i].x = 0.;
		frc[i].y = 0.;
		frc[i].z = 0.;
	}
}


static __global__ void MD_Iteration_Speed_Verlet_1(const int atom_numbers, const float half_dt, const float dt, const VECTOR *acc, VECTOR *vel, VECTOR *crd, VECTOR *frc)
{
	int i = blockDim.x*blockIdx.x + threadIdx.x;
	if (i < atom_numbers)
	{
		vel[i].x = vel[i].x + half_dt*acc[i].x;
		vel[i].y = vel[i].y + half_dt*acc[i].y;
		vel[i].z = vel[i].z + half_dt*acc[i].z;
		crd[i].x = crd[i].x + dt*vel[i].x;
		crd[i].y = crd[i].y + dt*vel[i].y;
		crd[i].z = crd[i].z + dt*vel[i].z;
		frc[i].x = 0.;
		frc[i].y = 0.;
		frc[i].z = 0.;
	}
}

static __global__ void MD_Iteration_Speed_Verlet_2(const int atom_numbers, const float half_dt, const float *inverse_mass, const VECTOR *frc, VECTOR *vel, VECTOR *acc)
{
	int i = blockDim.x*blockIdx.x + threadIdx.x;
	if (i < atom_numbers)
	{
		acc[i].x = inverse_mass[i] * frc[i].x;
		acc[i].y = inverse_mass[i] * frc[i].y;
		acc[i].z = inverse_mass[i] * frc[i].z;
		vel[i].x = vel[i].x + half_dt*acc[i].x;
		vel[i].y = vel[i].y + half_dt*acc[i].y;
		vel[i].z = vel[i].z + half_dt*acc[i].z;
	}
}

static __global__ void MD_Iteration_Speed_Verlet_2_With_Max_Velocity(const int atom_numbers, const float half_dt, const float *inverse_mass, const VECTOR *frc, VECTOR *vel, VECTOR *acc, const float max_velocity)
{
	int i = blockDim.x*blockIdx.x + threadIdx.x;
	if (i < atom_numbers)
	{
		VECTOR acc_i = inverse_mass[i] * frc[i];
		VECTOR vel_i = vel[i] + half_dt * acc_i;

		vel[i] = Make_Vector_Not_Exceed_Value(vel_i, max_velocity);
		acc[i] = acc_i;
	}
}


static __global__ void Get_Center_Of_Mass(const int residue_numbers, const int *start, const int *end,
	const VECTOR *crd, const float *atom_mass, const float *residue_mass_inverse, VECTOR *center_of_mass)
{
	for (int residue_i = blockDim.x*blockIdx.x + threadIdx.x; residue_i < residue_numbers; residue_i += gridDim.x * blockDim.x)
	{
		VECTOR com_lin = { 0.0f, 0.0f, 0.0f };
		for (int atom_i = start[residue_i]; atom_i < end[residue_i]; atom_i += 1)
		{
			com_lin = com_lin + atom_mass[atom_i] * crd[atom_i];
		}
		center_of_mass[residue_i] = residue_mass_inverse[residue_i] * com_lin;
	}
}

static __global__ void Map_Center_Of_Mass(const int residue_numbers, const int *start, const int *end,
	const float scaler, const VECTOR *center_of_mass, const VECTOR box_length, const VECTOR *no_wrap_crd, VECTOR *crd)
{
	VECTOR trans_vec;
	VECTOR com;
	for (int residue_i = blockDim.x*blockIdx.x + threadIdx.x; residue_i < residue_numbers; residue_i += gridDim.x * blockDim.x)
	{
		com = center_of_mass[residue_i];

		trans_vec.x = com.x - floorf(com.x / box_length.x) * box_length.x;
		trans_vec.y = com.y - floorf(com.y / box_length.y) * box_length.y;
		trans_vec.z = com.z - floorf(com.z / box_length.z) * box_length.z;
		trans_vec = scaler * trans_vec - com;

		for (int atom_i = start[residue_i] + threadIdx.y; atom_i < end[residue_i]; atom_i += blockDim.y)
		{
			crd[atom_i] = no_wrap_crd[atom_i] + trans_vec;
		}
	}
}

static __global__ void Add_Sum_List(int n, float *atom_virial, float *sum_virial)
{
	float temp = 0;
	for (int i = threadIdx.x; i < n; i = i + blockDim.x)
	{
		temp = temp + atom_virial[i];
	}
	atomicAdd(sum_virial, temp);
}

static __global__ void Calculate_Pressure_Cuda(const float V_inverse, const float *ek, const float *virial, float *pressure)
{
	pressure[0] = (ek[0] * 2 + virial[0]) * 0.33333333333333f * V_inverse;
}




static __global__ void MD_Temperature
(const int residue_numbers, const int *start, const int *end, float *ek,
const VECTOR *atom_vel, const float *atom_mass)
{
	int residue_i = blockDim.x*blockIdx.x + threadIdx.x;
	if (residue_i < residue_numbers)
	{
		VECTOR momentum = { 0., 0., 0. };
		float res_mass = 0.; //待提出，只需要初始时计算一遍
		int s = start[residue_i];
		int e = end[residue_i];
		float mass_lin;
		for (int atom_i = s; atom_i < e; atom_i = atom_i + 1)
		{
			mass_lin = atom_mass[atom_i];

			momentum.x = momentum.x + mass_lin*atom_vel[atom_i].x;
			momentum.y = momentum.y + mass_lin*atom_vel[atom_i].y;
			momentum.z = momentum.z + mass_lin*atom_vel[atom_i].z;
			res_mass = res_mass + mass_lin;
		}
		ek[residue_i] = 0.5*(momentum.x*momentum.x + momentum.y*momentum.y + momentum.z*momentum.z) / res_mass * 2. / 3. / CONSTANT_kB / residue_numbers;
	}
}

static __global__ void MD_Residue_Ek
(const int residue_numbers,const int *start,const int *end,float *ek,
const VECTOR *atom_vel,const float *atom_mass)
{
	int residue_i = blockDim.x*blockIdx.x + threadIdx.x;
	if (residue_i < residue_numbers)
	{
		VECTOR momentum = {0.,0.,0.};
		float res_mass = 0.; //待提出，只需要初始时计算一遍
		int s = start[residue_i];
		int e = end[residue_i];
		float mass_lin;
		for (int atom_i = s; atom_i < e; atom_i = atom_i + 1)
		{
			mass_lin = atom_mass[atom_i];

			momentum.x = momentum.x + mass_lin*atom_vel[atom_i].x;
			momentum.y = momentum.y + mass_lin*atom_vel[atom_i].y;
			momentum.z = momentum.z + mass_lin*atom_vel[atom_i].z;
			res_mass = res_mass + mass_lin;
		}
		ek[residue_i] = 0.5*(momentum.x*momentum.x + momentum.y*momentum.y + momentum.z*momentum.z) / res_mass;
	}
}

static __global__ void MD_Atom_Ek
(const int atom_numbers, float *ek, const VECTOR *atom_vel, const float *atom_mass)
{
	int atom_i = blockDim.x*blockIdx.x + threadIdx.x;
	if (atom_i < atom_numbers)
	{
		VECTOR v = atom_vel[atom_i];
		ek[atom_i] = 0.5 * v * v * atom_mass[atom_i];
	}
}

void MD_INFORMATION::system_information::Initial(CONTROLLER *controller, MD_INFORMATION *md_info)
{
	this->md_info = md_info;
	steps = 0;
	step_limit = 1000;
	if (controller[0].Command_Exist("step_limit"))
	{
		step_limit = atoi(controller[0].Command("step_limit"));
	}

	target_temperature = 300.0f;
	if (md_info->mode >= md_info->NVT && controller[0].Command_Exist("target_temperature"))
	{
		target_temperature = atof(controller[0].Command("target_temperature"));
	}

	target_pressure = 1;
	if (md_info->mode == md_info->NPT &&  controller[0].Command_Exist("target_pressure"))
		target_pressure = atof(controller[0].Command("target_pressure"));
	target_pressure *= CONSTANT_PRES_CONVERTION_INVERSE;

	controller->Step_Print_Initial("step", "%d");
	controller->Step_Print_Initial("time", "%.3lf");
	controller->Step_Print_Initial("temperature", "%.2f");
	controller->Step_Print_Initial("potential", "%.2f");
	Cuda_Malloc_Safely((void**)&this->d_virial, sizeof(float));
	Cuda_Malloc_Safely((void**)&this->d_pressure, sizeof(float));
	Cuda_Malloc_Safely((void**)&this->d_temperature, sizeof(float));
	Cuda_Malloc_Safely((void**)&this->d_potential, sizeof(float));
	Cuda_Malloc_Safely((void**)&this->d_sum_of_atom_ek, sizeof(float));
}

void MD_INFORMATION::non_bond_information::Initial(CONTROLLER *controller, MD_INFORMATION *md_info)
{
	if (controller[0].Command_Exist("skin"))
	{
		skin = atof(controller[0].Command("skin"));
	}
	else
	{
		skin = 2.0;
	}
	controller->printf("    skin set to %.2f Angstram\n", skin);

	if (controller[0].Command_Exist("cutoff"))
	{
		cutoff = atof(controller[0].Command("cutoff"));
	}
	else
	{
		cutoff = 10.0;
	}
	controller->printf("    cutoff set to %.2f Angstram\n", cutoff);
	/*===========================
	读取排除表相关信息
	============================*/
	if (controller[0].Command_Exist("exclude_in_file"))
	{
		FILE *fp = NULL;
		controller->printf("    Start reading excluded list:\n");
		Open_File_Safely(&fp, controller[0].Command("exclude_in_file"), "r");
		
		int atom_numbers = 0;
		fscanf(fp, "%d %d", &atom_numbers, &excluded_atom_numbers);
		if (md_info->atom_numbers > 0 && md_info->atom_numbers != atom_numbers)
		{
			controller->printf("        Error: atom_numbers is not equal: %d %d\n", md_info->atom_numbers, atom_numbers);
			getchar();
			exit(1);
		}
		else if (md_info->atom_numbers == 0)
		{
			md_info->atom_numbers = atom_numbers;
		}
		controller->printf("        excluded list total length is %d\n", excluded_atom_numbers);

		Cuda_Malloc_Safely((void**)&d_excluded_list_start, sizeof(int)*atom_numbers);
		Cuda_Malloc_Safely((void**)&d_excluded_numbers, sizeof(int)*atom_numbers);
		Cuda_Malloc_Safely((void**)&d_excluded_list, sizeof(int)*excluded_atom_numbers);

		Malloc_Safely((void**)&h_excluded_list_start, sizeof(int)*atom_numbers);
		Malloc_Safely((void**)&h_excluded_numbers, sizeof(int)*atom_numbers);
		Malloc_Safely((void**)&h_excluded_list, sizeof(int)*excluded_atom_numbers);
		int count = 0;
		for (int i = 0; i < atom_numbers; i++)
		{
			fscanf(fp, "%d", &h_excluded_numbers[i]);
			h_excluded_list_start[i] = count;
			for (int j = 0; j < h_excluded_numbers[i]; j++)
			{
				fscanf(fp, "%d", &h_excluded_list[count]);
				count++;
			}
		}
		cudaMemcpy(d_excluded_list_start, h_excluded_list_start, sizeof(int)*atom_numbers, cudaMemcpyHostToDevice);
		cudaMemcpy(d_excluded_numbers, h_excluded_numbers, sizeof(int)*atom_numbers, cudaMemcpyHostToDevice);
		cudaMemcpy(d_excluded_list, h_excluded_list, sizeof(int)*excluded_atom_numbers, cudaMemcpyHostToDevice);
		controller->printf("    End reading excluded list\n\n");
		fclose(fp);
	}
	else if (controller[0].Command_Exist("amber_parm7"))
	{
		/*===========================
		从parm中读取排除表相关信息
		============================*/
		FILE *parm = NULL;
		Open_File_Safely(&parm, controller[0].Command("amber_parm7"), "r");
		controller->printf("    Start reading excluded list from AMBER file:\n");
		while (true)
		{
			char temps[CHAR_LENGTH_MAX];
			char temp_first_str[CHAR_LENGTH_MAX];
			char temp_second_str[CHAR_LENGTH_MAX];
			if (!fgets(temps, CHAR_LENGTH_MAX, parm))
			{
				break;
			}
			if (sscanf(temps, "%s %s", temp_first_str, temp_second_str) != 2)
			{
				continue;
			}
			if (strcmp(temp_first_str, "%FLAG") == 0
				&& strcmp(temp_second_str, "POINTERS") == 0)
			{
				fgets(temps, CHAR_LENGTH_MAX, parm);

				int atom_numbers = 0;
				fscanf(parm, "%d\n", &atom_numbers);
				if (md_info->atom_numbers > 0 && md_info->atom_numbers != atom_numbers)
				{
					controller->printf("        Error: atom_numbers is not equal: %d %d\n", md_info->atom_numbers, atom_numbers);
					getchar();
					exit(1);
				}
				else if (md_info->atom_numbers == 0)
				{
					md_info->atom_numbers = atom_numbers;
				}
				Cuda_Malloc_Safely((void**)&d_excluded_list_start, sizeof(int)*atom_numbers);
				Cuda_Malloc_Safely((void**)&d_excluded_numbers, sizeof(int)*atom_numbers);

				Malloc_Safely((void**)&h_excluded_list_start, sizeof(int)*atom_numbers);
				Malloc_Safely((void**)&h_excluded_numbers, sizeof(int)*atom_numbers);
				for (int i = 0; i < 9; i = i + 1)
				{
					fscanf(parm, "%d\n", &excluded_atom_numbers);
				}
				fscanf(parm, "%d\n", &excluded_atom_numbers);
				controller->printf("        excluded list total length is %d\n", excluded_atom_numbers);

				Cuda_Malloc_Safely((void**)&d_excluded_list, sizeof(int)*excluded_atom_numbers);
				Malloc_Safely((void**)&h_excluded_list, sizeof(int)*excluded_atom_numbers);
			}

			//read atom_excluded_number for every atom
			if (strcmp(temp_first_str, "%FLAG") == 0
				&& strcmp(temp_second_str, "NUMBER_EXCLUDED_ATOMS") == 0)
			{
				fgets(temps, CHAR_LENGTH_MAX, parm);
				for (int i = 0; i<md_info->atom_numbers; i = i + 1)
				{
					fscanf(parm, "%d\n", &h_excluded_numbers[i]);
				}
			}
			//read every atom's excluded atom list
			if (strcmp(temp_first_str, "%FLAG") == 0
				&& strcmp(temp_second_str, "EXCLUDED_ATOMS_LIST") == 0)
			{
				int count = 0;
				int none_count = 0;
				int lin = 0;
				fgets(temps, CHAR_LENGTH_MAX, parm);
				for (int i = 0; i<md_info->atom_numbers; i = i + 1)
				{
					h_excluded_list_start[i] = count;
					for (int j = 0; j<h_excluded_numbers[i]; j = j + 1)
					{
						fscanf(parm, "%d\n", &lin);
						if (lin == 0)
						{
							h_excluded_numbers[i] = 0;
							break;
						}
						else
						{
							h_excluded_list[count] = lin - 1;
							count = count + 1;
						}
					}
					if (h_excluded_numbers[i] > 0)
						thrust::sort(&h_excluded_list[h_excluded_list_start[i]], &h_excluded_list[h_excluded_list_start[i]] + h_excluded_numbers[i]);
				}
			}
		}

		cudaMemcpy(d_excluded_list_start, h_excluded_list_start, sizeof(int)*md_info->atom_numbers, cudaMemcpyHostToDevice);
		cudaMemcpy(d_excluded_numbers, h_excluded_numbers, sizeof(int)*md_info->atom_numbers, cudaMemcpyHostToDevice);
		cudaMemcpy(d_excluded_list, h_excluded_list, sizeof(int)*excluded_atom_numbers, cudaMemcpyHostToDevice);
		controller->printf("    End reading excluded list from AMBER file\n\n");
		fclose(parm);
	}
	else
	{
		int atom_numbers = md_info->atom_numbers;
		excluded_atom_numbers = 0;
		controller->printf("    Set all atom exclude no atoms as default\n"); 

		Cuda_Malloc_Safely((void**)&d_excluded_list_start, sizeof(int)*atom_numbers);
		Cuda_Malloc_Safely((void**)&d_excluded_numbers, sizeof(int)*atom_numbers);
		Cuda_Malloc_Safely((void**)&d_excluded_list, sizeof(int)*excluded_atom_numbers);

		Malloc_Safely((void**)&h_excluded_list_start, sizeof(int)*atom_numbers);
		Malloc_Safely((void**)&h_excluded_numbers, sizeof(int)*atom_numbers);
		Malloc_Safely((void**)&h_excluded_list, sizeof(int)*excluded_atom_numbers);


		int count = 0;
		for (int i = 0; i < atom_numbers; i++)
		{
			h_excluded_numbers[i] = 0;
			h_excluded_list_start[i] = count;
			for (int j = 0; j < h_excluded_numbers[i]; j++)
			{
				h_excluded_list[count] = 0;
				count++;
			}
		}
		cudaMemcpy(d_excluded_list_start, h_excluded_list_start, sizeof(int)*atom_numbers, cudaMemcpyHostToDevice);
		cudaMemcpy(d_excluded_numbers, h_excluded_numbers, sizeof(int)*atom_numbers, cudaMemcpyHostToDevice);
		cudaMemcpy(d_excluded_list, h_excluded_list, sizeof(int)*excluded_atom_numbers, cudaMemcpyHostToDevice);
	}

}

void MD_INFORMATION::periodic_box_condition_information::Initial(CONTROLLER *controller, VECTOR box_length)
{
	crd_to_uint_crd_cof = CONSTANT_UINT_MAX_FLOAT / box_length;
	quarter_crd_to_uint_crd_cof = 0.25 * crd_to_uint_crd_cof;
	uint_dr_to_dr_cof = 1.0f / crd_to_uint_crd_cof;
}

void MD_INFORMATION::Read_Mode(CONTROLLER *controller)
{
	if (controller[0].Command_Exist("mode"))
	{
		if (is_str_equal(controller[0].Command("mode"), "NVT"))
		{
			controller->printf("    Mode set to NVT\n");
			mode = 1;
		}
		else if (is_str_equal(controller[0].Command("mode"), "NPT"))
		{
			controller->printf("    Mode set to NPT\n");
			mode = 2;
		}
		else if (is_str_equal(controller[0].Command("mode"), "Minimization"))
		{
			controller->printf("    Mode set to Energy Minimization\n");
			mode = -1;
		}
		else if (is_str_equal(controller[0].Command("mode"), "NVE"))
		{
			controller->printf("    Mode set to NVE\n");
			mode = 0;
		}
		else
		{
			controller->printf("    Warning: Mode '%s' not match. Set to NVE as default\n", controller[0].Command("mode"));
			mode = 0;
		}
	}
	else
	{
		controller->printf("    Mode set to NVE as default\n");
		mode = 0;
	}
}

void MD_INFORMATION::Read_dt(CONTROLLER *controller)
{
	if (controller[0].Command_Exist("dt"))
	{
		controller->printf("    dt set to %f ps\n", atof(controller[0].Command("dt")));
		dt = atof(controller[0].Command("dt")) * CONSTANT_TIME_CONVERTION;
		sscanf(controller[0].Command("dt"), "%lf", &sys.dt_in_ps);
	}
	else
	{
		dt = 0.001 * CONSTANT_TIME_CONVERTION;
		sys.dt_in_ps = 0.001;
		controller->printf("    dt set to %f ps\n", 0.001);
	}
}

void MD_INFORMATION::trajectory_output::Initial(CONTROLLER *controller, MD_INFORMATION *md_info)
{
	this->md_info = md_info;
	current_crd_synchronized_step = 0;
	is_molecule_map_output = 0;
	if (controller[0].Command_Exist("molecule_map_output"))
	{
		is_molecule_map_output = atoi(controller[0].Command("molecule_map_output"));
	}
	write_trajectory_interval = 1000;
	if (controller[0].Command_Exist("write_information_interval"))
	{
		write_trajectory_interval = atoi(controller[0].Command("write_information_interval"));
	}
	write_restart_file_interval = write_trajectory_interval;
	if (controller[0].Command_Exist("write_restart_file_interval"))
	{
		write_restart_file_interval = atoi(controller[0].Command("write_restart_file_interval"));
	}
	if (controller->Command_Exist(TRAJ_COMMAND))
	{
		Open_File_Safely(&crd_traj, controller->Command(TRAJ_COMMAND), "wb");
	}
	else
	{
		Open_File_Safely(&crd_traj, TRAJ_DEFAULT_FILENAME, "wb");
	}
	if (controller->Command_Exist(BOX_TRAJ_COMMAND))
	{
		Open_File_Safely(&box_traj, controller->Command(BOX_TRAJ_COMMAND), "w");
	}
	else
	{
		Open_File_Safely(&box_traj, BOX_TRAJ_DEFAULT_FILENAME, "w");
	}
	if (controller->Command_Exist(RESTART_COMMAND))
	{
		strcpy(restart_name, controller->Command(RESTART_COMMAND));
	}
	else
	{
		strcpy(restart_name, RESTART_DEFAULT_FILENAME);
	}
}

void MD_INFORMATION::NVE_iteration::Initial(CONTROLLER *controller, MD_INFORMATION *md_info)
{
	this->md_info = md_info;
	max_velocity = -1;
	if (controller[0].Command_Exist("nve_velocity_max"))
	{
		max_velocity = atof(controller[0].Command("nve_velocity_max"));
	}
}
void MD_INFORMATION::residue_information::Read_AMBER_Parm7(const char *file_name, CONTROLLER controller)
{
	FILE *parm = NULL;
	Open_File_Safely(&parm, file_name, "r");
	controller.printf("    Start reading residue informataion from AMBER parm7:\n");
	
	while (true)
	{
		char temps[CHAR_LENGTH_MAX];
		char temp_first_str[CHAR_LENGTH_MAX];
		char temp_second_str[CHAR_LENGTH_MAX];
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

			int atom_numbers = 0;
			fscanf(parm, "%d", &atom_numbers);
			if (md_info->atom_numbers > 0 && md_info->atom_numbers != atom_numbers)
			{
				controller.printf("        Error: atom_numbers is not equal: %d %d\n", md_info->atom_numbers, atom_numbers);
				getchar();
				exit(1);
			}
			else if (md_info->atom_numbers == 0)
			{
				md_info->atom_numbers = atom_numbers;
			}
			for (int i = 0; i < 10; i = i + 1)
			{
				int lin;
				fscanf(parm, "%d\n", &lin);
			}
			fscanf(parm, "%d\n", &this->residue_numbers);//NRES
			controller.printf("        residue_numbers is %d\n", this->residue_numbers);

			Malloc_Safely((void**)&h_mass, sizeof(float)*this->residue_numbers);
			Malloc_Safely((void**)&h_mass_inverse, sizeof(float)*this->residue_numbers);
			Malloc_Safely((void**)&h_res_start, sizeof(int)*this->residue_numbers);
			Malloc_Safely((void**)&h_res_end, sizeof(int)*this->residue_numbers);
			Malloc_Safely((void**)&h_momentum, sizeof(float)*this->residue_numbers);
			Malloc_Safely((void**)&h_center_of_mass, sizeof(VECTOR)*this->residue_numbers);
			Malloc_Safely((void**)&h_sigma_of_res_ek, sizeof(float));

			Cuda_Malloc_Safely((void**)&d_mass, sizeof(float)*this->residue_numbers);
			Cuda_Malloc_Safely((void**)&d_mass_inverse, sizeof(float)*this->residue_numbers);
			Cuda_Malloc_Safely((void**)&d_res_start, sizeof(int)*this->residue_numbers);
			Cuda_Malloc_Safely((void**)&d_res_end, sizeof(int)*this->residue_numbers);
			Cuda_Malloc_Safely((void**)&d_momentum, sizeof(float)*this->residue_numbers);
			Cuda_Malloc_Safely((void**)&d_center_of_mass, sizeof(VECTOR)*this->residue_numbers);
			Cuda_Malloc_Safely((void**)&res_ek_energy, sizeof(float)*this->residue_numbers); 
			Cuda_Malloc_Safely((void**)&sigma_of_res_ek, sizeof(float));
		}//FLAG POINTERS

		//residue range read
		if (strcmp(temp_first_str, "%FLAG") == 0
			&& strcmp(temp_second_str, "RESIDUE_POINTER") == 0)
		{
			fgets(temps, CHAR_LENGTH_MAX, parm);
			//注意读进来的数的编号要减1
			int *lin_serial;
			Malloc_Safely((void**)&lin_serial, sizeof(int)* this->residue_numbers);
			for (int i = 0; i<this->residue_numbers; i = i + 1)
			{
				fscanf(parm, "%d\n", &lin_serial[i]);
			}
			for (int i = 0; i<this->residue_numbers - 1; i = i + 1)
			{
				h_res_start[i] = lin_serial[i] - 1;
				h_res_end[i] = lin_serial[i + 1] - 1;
			}
			h_res_start[this->residue_numbers - 1] = lin_serial[this->residue_numbers - 1] - 1;
			h_res_end[this->residue_numbers - 1] = md_info->atom_numbers + 1 - 1;

			free(lin_serial);
		}
	}//while cycle

	cudaMemcpy(this->d_res_start, h_res_start, sizeof(int)*this->residue_numbers, cudaMemcpyHostToDevice);
	cudaMemcpy(this->d_res_end, h_res_end, sizeof(int)*this->residue_numbers, cudaMemcpyHostToDevice);

	controller.printf("    End reading residue informataion from AMBER parm7\n");

	fclose(parm);
}

void MD_INFORMATION::residue_information::Initial(CONTROLLER *controller, MD_INFORMATION *md_info)
{
	this->md_info = md_info;
	if (!(controller[0].Command_Exist("residue_in_file")))
	{
		if (controller[0].Command_Exist("amber_parm7"))
		{
			Read_AMBER_Parm7(controller[0].Command("amber_parm7"), controller[0]);
			is_initialzed = 1;
		}
		//对于没有residue输入的模拟，默认每个粒子作为一个residue
		else
		{
			residue_numbers = md_info->atom_numbers;
			controller->printf("    Set default residue list:\n");
			controller->printf("        residue_numbers is %d\n", residue_numbers);
			Malloc_Safely((void**)&h_mass, sizeof(float)*this->residue_numbers);
			Malloc_Safely((void**)&h_mass_inverse, sizeof(float)*this->residue_numbers);
			Malloc_Safely((void**)&h_res_start, sizeof(int)*this->residue_numbers);
			Malloc_Safely((void**)&h_res_end, sizeof(int)*this->residue_numbers);
			Malloc_Safely((void**)&h_momentum, sizeof(float)*this->residue_numbers);
			Malloc_Safely((void**)&h_center_of_mass, sizeof(VECTOR)*this->residue_numbers);
			Malloc_Safely((void**)&h_sigma_of_res_ek, sizeof(float));

			Cuda_Malloc_Safely((void**)&d_mass, sizeof(float)*this->residue_numbers);
			Cuda_Malloc_Safely((void**)&d_mass_inverse, sizeof(float)*this->residue_numbers);
			Cuda_Malloc_Safely((void**)&d_res_start, sizeof(int)*this->residue_numbers);
			Cuda_Malloc_Safely((void**)&d_res_end, sizeof(int)*this->residue_numbers);
			Cuda_Malloc_Safely((void**)&d_momentum, sizeof(float)*this->residue_numbers);
			Cuda_Malloc_Safely((void**)&d_center_of_mass, sizeof(VECTOR)*this->residue_numbers);
			Cuda_Malloc_Safely((void**)&res_ek_energy, sizeof(float)*this->residue_numbers);
			Cuda_Malloc_Safely((void**)&sigma_of_res_ek, sizeof(float));
			int count = 0;
			int temp=1;//每个粒子作为一个residue
			for (int i = 0; i < residue_numbers; i++)
			{
				h_res_start[i] = count;
				count += temp;
				h_res_end[i] = count;
			}
			cudaMemcpy(d_res_start, h_res_start, sizeof(int)* residue_numbers, cudaMemcpyHostToDevice);
			cudaMemcpy(d_res_end, h_res_end, sizeof(int)* residue_numbers, cudaMemcpyHostToDevice);
			controller->printf("    End reading residue list\n\n");
			is_initialzed = 1;
		}
	}
	else
	{
		FILE *fp = NULL;
		controller->printf("    Start reading residue list:\n");
		Open_File_Safely(&fp, controller[0].Command("residue_in_file"), "r");
		int atom_numbers = 0;
		fscanf(fp, "%d %d", &atom_numbers, &residue_numbers);
		if (md_info->atom_numbers > 0 && md_info->atom_numbers != atom_numbers)
		{
			controller->printf("        Error: atom_numbers is not equal: %d %d\n", md_info->atom_numbers, atom_numbers);
			getchar();
			exit(1);
		}
		else if (md_info->atom_numbers == 0)
		{
			md_info->atom_numbers = atom_numbers;
		}
		controller->printf("        residue_numbers is %d\n", residue_numbers);
		Malloc_Safely((void**)&h_mass, sizeof(float)*this->residue_numbers);
		Malloc_Safely((void**)&h_mass_inverse, sizeof(float)*this->residue_numbers);
		Malloc_Safely((void**)&h_res_start, sizeof(int)*this->residue_numbers);
		Malloc_Safely((void**)&h_res_end, sizeof(int)*this->residue_numbers);
		Malloc_Safely((void**)&h_momentum, sizeof(float)*this->residue_numbers);
		Malloc_Safely((void**)&h_center_of_mass, sizeof(VECTOR)*this->residue_numbers);
		Malloc_Safely((void**)&h_sigma_of_res_ek, sizeof(float));

		Cuda_Malloc_Safely((void**)&d_mass, sizeof(float)*this->residue_numbers);
		Cuda_Malloc_Safely((void**)&d_mass_inverse, sizeof(float)*this->residue_numbers);
		Cuda_Malloc_Safely((void**)&d_res_start, sizeof(int)*this->residue_numbers);
		Cuda_Malloc_Safely((void**)&d_res_end, sizeof(int)*this->residue_numbers);
		Cuda_Malloc_Safely((void**)&d_momentum, sizeof(float)*this->residue_numbers);
		Cuda_Malloc_Safely((void**)&d_center_of_mass, sizeof(VECTOR)*this->residue_numbers);
		Cuda_Malloc_Safely((void**)&res_ek_energy, sizeof(float)*this->residue_numbers);
		Cuda_Malloc_Safely((void**)&sigma_of_res_ek, sizeof(float));

		int count = 0;
		int temp;
		for (int i = 0; i < residue_numbers; i++)
		{
			h_res_start[i] = count;
			fscanf(fp, "%d", &temp);
			count += temp;
			h_res_end[i] = count;
		}
		cudaMemcpy(d_res_start, h_res_start, sizeof(int)* residue_numbers, cudaMemcpyHostToDevice);
		cudaMemcpy(d_res_end, h_res_end, sizeof(int)* residue_numbers, cudaMemcpyHostToDevice);
		controller->printf("    End reading residue list\n\n");
		fclose(fp);
		is_initialzed = 1;
	}
	if (is_initialzed)
	{
		if (md_info->h_mass != NULL)
		{
			for (int i = 0; i < residue_numbers; i++)
			{
				float temp_mass = 0;
				for (int j = h_res_start[i]; j < h_res_end[i]; j++)
				{
					temp_mass += md_info->h_mass[j];
				}
				this->h_mass[i] = temp_mass;
				if (temp_mass == 0)
					this->h_mass_inverse[i] = 0;
				else
					this->h_mass_inverse[i] = 1.0 / temp_mass;
			}
			cudaMemcpy(d_mass_inverse, h_mass_inverse, sizeof(float)* residue_numbers, cudaMemcpyHostToDevice);
			cudaMemcpy(d_mass, h_mass, sizeof(float)* residue_numbers, cudaMemcpyHostToDevice);
		}
		else
		{
			controller->printf("    Error: atom mass should be initialized before residue mass\n");
			getchar();
			exit(1);
		}
	}
}

void MD_INFORMATION::Read_Coordinate_And_Velocity(CONTROLLER *controller)
{
	sys.start_time = 0.0;
	if (controller[0].Command_Exist("coordinate_in_file"))
	{
		Read_Coordinate_In_File(controller[0].Command("coordinate_in_file"), controller[0]);
		if (controller[0].Command_Exist("velocity_in_file"))
		{
			FILE *fp = NULL;
			controller->printf("    Start reading velocity_in_file:\n");
			Open_File_Safely(&fp, controller[0].Command("velocity_in_file"), "r");
			
			int atom_numbers = 0;
			char lin[CHAR_LENGTH_MAX];
			fgets(lin, CHAR_LENGTH_MAX, fp);
			int scanf_ret = sscanf(lin, "%d", &atom_numbers);
			if (this->atom_numbers > 0 && this->atom_numbers != atom_numbers)
			{
				controller->printf("        Error: atom_numbers is not equal: %d %d\n", this->atom_numbers, atom_numbers);
				getchar();
				exit(1);
			}
			Malloc_Safely((void**)&velocity, sizeof(VECTOR)*this->atom_numbers);
			Cuda_Malloc_Safely((void**)&vel, sizeof(VECTOR)*this->atom_numbers);
			for (int i = 0; i < atom_numbers; i++)
			{
				fscanf(fp, "%f %f %f", &velocity[i].x, &velocity[i].y, &velocity[i].z);
			}
			cudaMemcpy(vel, velocity, sizeof(VECTOR)* atom_numbers, cudaMemcpyHostToDevice);
			controller->printf("    End reading velocity_in_file\n\n");
			fclose(fp);
		}
		else
		{
			controller->printf("    Velocity is set to zero as default\n");
			Malloc_Safely((void**)&velocity, sizeof(VECTOR)*this->atom_numbers);
			Cuda_Malloc_Safely((void**)&vel, sizeof(VECTOR)*this->atom_numbers);
			for (int i = 0; i < atom_numbers; i++)
			{
				velocity[i].x = 0;
				velocity[i].y = 0;
				velocity[i].z = 0;
			}
			cudaMemcpy(vel, velocity, sizeof(VECTOR)* atom_numbers, cudaMemcpyHostToDevice);
		}
	}
	else if (controller[0].Command_Exist("amber_rst7"))
	{
		output.amber_irest = 1;
		if (controller[0].Command_Exist("amber_irest"))
			output.amber_irest = atoi(controller[0].Command("amber_irest"));
		Read_Rst7(controller[0].Command("amber_rst7"), output.amber_irest, controller[0]);
	}
	else
	{
		printf("MD basic information needed. Specify the coordinate in file.\n");
		getchar();
		exit(1);
	}
}

void MD_INFORMATION::Read_Mass(CONTROLLER *controller)
{
	if (controller[0].Command_Exist("mass_in_file"))
	{
		FILE *fp = NULL;
		controller->printf("    Start reading mass:\n");
		Open_File_Safely(&fp, controller[0].Command("mass_in_file"), "r");
		int atom_numbers = 0;
		char lin[CHAR_LENGTH_MAX];
		fgets(lin, CHAR_LENGTH_MAX, fp);
		int scanf_ret = sscanf(lin, "%d", &atom_numbers);
		if (this->atom_numbers > 0 && this->atom_numbers != atom_numbers)
		{
			controller->printf("        Error: atom_numbers is not equal: %d %d\n", this->atom_numbers, atom_numbers);
			getchar();
			exit(1);
		}
		else if (this->atom_numbers == 0)
		{
			this->atom_numbers = atom_numbers;
		}
		Malloc_Safely((void**)&h_mass, sizeof(float)* atom_numbers);
		Malloc_Safely((void**)&h_mass_inverse, sizeof(float)* atom_numbers);
		Cuda_Malloc_Safely((void**)&d_mass, sizeof(float)* atom_numbers);
		Cuda_Malloc_Safely((void**)&d_mass_inverse, sizeof(float)* atom_numbers);
		sys.total_mass = 0;
		for (int i = 0; i < atom_numbers; i++)
		{
			fscanf(fp, "%f", &h_mass[i]);
			sys.total_mass += h_mass[i];
			if (h_mass[i] == 0)
				h_mass_inverse[i] = 0;
			else
				h_mass_inverse[i] = 1.0 / h_mass[i];

		}
		controller->printf("    End reading mass\n\n");
		fclose(fp);
	}
	else if (controller[0].Command_Exist("amber_parm7"))
	{
		FILE *parm = NULL;
		Open_File_Safely(&parm, controller[0].Command("amber_parm7"), "r");
		controller[0].printf("    Start reading mass from AMBER parm7:\n");
		while (true)
		{
			char temps[CHAR_LENGTH_MAX];
			char temp_first_str[CHAR_LENGTH_MAX];
			char temp_second_str[CHAR_LENGTH_MAX];
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

				int atom_numbers = 0;
				fscanf(parm, "%d", &atom_numbers);
				if (this->atom_numbers > 0 && this->atom_numbers != atom_numbers)
				{
					controller->printf("        Error: atom_numbers is not equal: %d %d\n", this->atom_numbers, atom_numbers);
					getchar();
					exit(1);
				}
				else if (this->atom_numbers == 0)
				{
					this->atom_numbers = atom_numbers;
				}
				Malloc_Safely((void**)&h_mass, sizeof(float)* atom_numbers);
				Malloc_Safely((void**)&h_mass_inverse, sizeof(float)* atom_numbers);
				Cuda_Malloc_Safely((void**)&d_mass, sizeof(float)* atom_numbers);
				Cuda_Malloc_Safely((void**)&d_mass_inverse, sizeof(float)* atom_numbers);
			}
			if (strcmp(temp_first_str, "%FLAG") == 0
				&& strcmp(temp_second_str, "MASS") == 0)
			{
				fgets(temps, CHAR_LENGTH_MAX, parm);
				double lin;
				sys.total_mass = 0;
				for (int i = 0; i < this->atom_numbers; i = i + 1)
				{
					fscanf(parm, "%lf\n", &lin);
					this->h_mass[i] = (float)lin;
					if (h_mass[i] == 0)
						h_mass_inverse[i] = 0;
					else
						h_mass_inverse[i] = 1.0f / h_mass[i];
					sys.total_mass += h_mass[i];
				}
			}
		}
		controller[0].printf("    End reading mass from AMBER parm7\n\n");
		fclose(parm);
	}
	else if (atom_numbers > 0)
	{
		controller[0].printf("    mass is set to 20 as default\n");
		sys.total_mass = 0;
		Malloc_Safely((void**)&h_mass, sizeof(float)* atom_numbers);
		Malloc_Safely((void**)&h_mass_inverse, sizeof(float)* atom_numbers);
		Cuda_Malloc_Safely((void**)&d_mass, sizeof(float)* atom_numbers);
		Cuda_Malloc_Safely((void**)&d_mass_inverse, sizeof(float)* atom_numbers);
		for (int i = 0; i < atom_numbers; i++)
		{
			h_mass[i] = 20;
			h_mass_inverse[i] = 1.0 / h_mass[i];
			sys.total_mass += h_mass[i];
		}
	}
	else
	{
		controller[0].printf("    Error: failed to initialze mass, because no atom_numbers found\n");
		getchar();
		exit(1);
	}
	if (atom_numbers > 0)
	{
		cudaMemcpy(d_mass, h_mass, sizeof(float)* atom_numbers, cudaMemcpyHostToDevice);
		cudaMemcpy(d_mass_inverse, h_mass_inverse, sizeof(float)* atom_numbers, cudaMemcpyHostToDevice);
	}
}

void MD_INFORMATION::Read_Charge(CONTROLLER *controller)
{
	if (controller[0].Command_Exist("charge_in_file"))
	{
		FILE *fp = NULL;
		controller->printf("    Start reading charge:\n");
		Open_File_Safely(&fp, controller[0].Command("charge_in_file"), "r");
		int atom_numbers = 0;
		char lin[CHAR_LENGTH_MAX];
		fgets(lin, CHAR_LENGTH_MAX, fp);
		int scanf_ret = sscanf(lin, "%d", &atom_numbers);
		if (this->atom_numbers > 0 && this->atom_numbers != atom_numbers)
		{
			controller->printf("        Error: atom_numbers is not equal: %d %d\n", this->atom_numbers, atom_numbers);
			getchar();
			exit(1);
		}
		else if (this->atom_numbers == 0)
		{
			this->atom_numbers = atom_numbers;
		}
		Malloc_Safely((void**)&h_charge, sizeof(float)* atom_numbers);
		Cuda_Malloc_Safely((void**)&d_charge, sizeof(float)* atom_numbers);
		for (int i = 0; i < atom_numbers; i++)
		{
			fscanf(fp, "%f", &h_charge[i]);
		}
		controller->printf("    End reading charge\n\n");
		fclose(fp);
	}
	else if (controller[0].Command_Exist("amber_parm7"))
	{
		FILE *parm = NULL;
		Open_File_Safely(&parm, controller[0].Command("amber_parm7"), "r");
		controller[0].printf("    Start reading charge from AMBER parm7:\n");
		while (true)
		{
			char temps[CHAR_LENGTH_MAX];
			char temp_first_str[CHAR_LENGTH_MAX];
			char temp_second_str[CHAR_LENGTH_MAX];
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

				int atom_numbers = 0;
				fscanf(parm, "%d", &atom_numbers);
				if (this->atom_numbers > 0 && this->atom_numbers != atom_numbers)
				{
					controller->printf("        Error: atom_numbers is not equal: %d %d\n", this->atom_numbers, atom_numbers);
					getchar();
					exit(1);
				}
				else if (this->atom_numbers == 0)
				{
					this->atom_numbers = atom_numbers;
				}
				Malloc_Safely((void**)&h_charge, sizeof(float)* atom_numbers);
				Cuda_Malloc_Safely((void**)&d_charge, sizeof(float)* atom_numbers);
			}
			if (strcmp(temp_first_str, "%FLAG") == 0
				&& strcmp(temp_second_str, "CHARGE") == 0)
			{
				fgets(temps, CHAR_LENGTH_MAX, parm);
				double lin;
				for (int i = 0; i < this->atom_numbers; i = i + 1)
				{
					fscanf(parm, "%f", &h_charge[i]);
				}
			}
		}
		controller[0].printf("    End reading charge from AMBER parm7\n\n");
		fclose(parm);
	}
	else if (atom_numbers > 0)
	{
		controller[0].printf("    charge is set to 0 as default\n");
		Malloc_Safely((void**)&h_charge, sizeof(float)* atom_numbers);
		Cuda_Malloc_Safely((void**)&d_charge, sizeof(float)* atom_numbers);
		for (int i = 0; i < atom_numbers; i++)
		{
			h_charge[i] = 0;
		}
	}
	else
	{
		controller[0].printf("    Error: failed to initialze charge, because no atom_numbers found\n");
		getchar();
		exit(1);
	}
	if (atom_numbers > 0)
	{
		cudaMemcpy(d_charge, h_charge, sizeof(float)* atom_numbers, cudaMemcpyHostToDevice);
	}
}

//MD_INFORMATION成员函数
void MD_INFORMATION::Initial(CONTROLLER *controller)
{
	controller->printf("START INITIALZING MD CORE:\n");
	atom_numbers = 0; //初始化，使得能够进行所有原子数目是否相等的判断
	
	strcpy(md_name, controller[0].Command("md_name"));
	Read_Mode(controller);
	Read_dt(controller);

	Read_Coordinate_And_Velocity(controller);
	
	
	Read_Mass(controller);
	Read_Charge(controller);
    
	sys.Initial(controller, this);  //!需要先初始化坐标和速度
	nb.Initial(controller, this);
	
	output.Initial(controller, this);

	nve.Initial(controller, this);
	
	res.Initial(controller, this);

	pbc.Initial(controller, sys.box_length);
	
	Atom_Information_Initial();

	is_initialized = 1;
	controller->printf("    structure last modify date is %d\n", last_modify_date);
	controller->printf("END INITIALZING MD CORE\n\n");
}

void MD_INFORMATION::Atom_Information_Initial()
{
	Malloc_Safely((void**)&this->force, sizeof(VECTOR)*this->atom_numbers);
	Malloc_Safely((void**)&this->h_atom_energy, sizeof(float)* atom_numbers);
	Malloc_Safely((void**)&this->h_atom_virial, sizeof(double)* atom_numbers);
	Cuda_Malloc_Safely((void**)&this->acc, sizeof(VECTOR)*this->atom_numbers);
	Cuda_Malloc_Safely((void**)&this->frc, sizeof(VECTOR)*this->atom_numbers);
	Cuda_Malloc_Safely((void **)&this->uint_crd, sizeof(UNSIGNED_INT_VECTOR)*this->atom_numbers);
	Cuda_Malloc_Safely((void**)&this->d_atom_energy, sizeof(float)* atom_numbers);
	Cuda_Malloc_Safely((void**)&this->d_atom_virial, sizeof(float)* atom_numbers);
	Cuda_Malloc_Safely((void**)&this->d_atom_ek, sizeof(float)* atom_numbers);
	Reset_List << <ceilf((float)3.*this->atom_numbers / 32), 32 >> >
		(3 * this->atom_numbers, (float*)this->acc, 0.);
	Reset_List << <ceilf((float)3.*this->atom_numbers / 32), 32 >> >
		(3 * this->atom_numbers, (float*)this->frc, 0.);
	sys.freedom = 3 * atom_numbers; //最大自由度，后面减
}

void MD_INFORMATION::Read_Coordinate_In_File(const char* file_name, CONTROLLER controller)
{
	FILE *fp =NULL;
	controller.printf("    Start reading coordinate_in_file:\n");
	Open_File_Safely(&fp, file_name, "r");
	char lin[CHAR_LENGTH_MAX];
	fgets(lin, CHAR_LENGTH_MAX, fp);
	int atom_numbers = 0;
	int scanf_ret = sscanf(lin, "%d %lf", &atom_numbers, &sys.start_time);
	if (this->atom_numbers > 0 && this->atom_numbers != atom_numbers)
	{
		controller.printf("        Error: atom_numbers is not equal: %d %d\n", this->atom_numbers, atom_numbers);
		getchar();
		exit(1);
	}
	else if (this->atom_numbers == 0)
	{
		this->atom_numbers = atom_numbers;
	}
	if (scanf_ret == 0)
	{
		controller.printf("        Error: Atom_numbers not found.\n");
		getchar();
		exit(1);
	}
	else  if (scanf_ret == 1)
	{
		sys.start_time = 0;
	}

	controller.printf("        atom_numbers is %d\n", this->atom_numbers);
	controller.printf("        system start_time is %lf\n", this->sys.start_time);
	Malloc_Safely((void**)&coordinate, sizeof(VECTOR)*this->atom_numbers);
	Cuda_Malloc_Safely((void**)&crd, sizeof(VECTOR)*this->atom_numbers);

	for (int i = 0; i < atom_numbers; i++)
	{
		fscanf(fp, "%f %f %f", &coordinate[i].x, &coordinate[i].y, &coordinate[i].z);
	}
	fscanf(fp, "%f %f %f", &sys.box_length.x, &sys.box_length.y, &sys.box_length.z);
	controller.printf("        box_length is\n            x: %f\n            y: %f\n            z: %f\n", sys.box_length.x, sys.box_length.y, sys.box_length.z);
	cudaMemcpy(crd, coordinate, sizeof(VECTOR)* atom_numbers, cudaMemcpyHostToDevice);
	controller.printf("    End reading coordinate_in_file\n\n");
	fclose(fp);
}
void MD_INFORMATION::Read_Rst7(const char* file_name, int irest, CONTROLLER controller)
{
	FILE *fin = NULL;
	Open_File_Safely(&fin, file_name, "r");
	controller.printf("    Start reading AMBER rst7:\n");
	char lin[CHAR_LENGTH_MAX];
	int atom_numbers = 0;
	fgets(lin, CHAR_LENGTH_MAX, fin);
	fgets(lin, CHAR_LENGTH_MAX, fin);
	int has_vel = 0;
	int scanf_ret = sscanf(lin, "%d %lf", &atom_numbers, &sys.start_time);
	if (this->atom_numbers > 0 && this->atom_numbers != atom_numbers)
	{
		controller.printf("        Error: atom_numbers is not equal: %d %d\n", this->atom_numbers, atom_numbers);
		getchar();
		exit(1);
	}
	else if (this->atom_numbers == 0)
	{
		this->atom_numbers = atom_numbers;
	}
	if (scanf_ret == 0)
	{
		controller.printf("        Error: Atom_numbers not found.\n");
		getchar();
		exit(1);
	}
	else  if (scanf_ret == 2)
	{
		has_vel = 1;
	}
	else
	{
		sys.start_time = 0;
	}

	Malloc_Safely((void**)&coordinate, sizeof(VECTOR)*this->atom_numbers);
	Cuda_Malloc_Safely((void**)&crd, sizeof(VECTOR)*this->atom_numbers);
	Malloc_Safely((void**)&velocity, sizeof(VECTOR)*this->atom_numbers);
	Cuda_Malloc_Safely((void**)&vel, sizeof(VECTOR)*this->atom_numbers);

	controller.printf("        atom_numbers is %d\n", this->atom_numbers);
	controller.printf("        system start time is %lf\n", this->sys.start_time);

	if (has_vel == 0 || irest == 0)
	{
		controller.printf("        All velocity will be set to 0\n");
	}


	for (int i = 0; i < this->atom_numbers; i = i + 1)
	{
		fscanf(fin, "%f %f %f",
			&this->coordinate[i].x,
			&this->coordinate[i].y,
			&this->coordinate[i].z);
	}
	if (has_vel)
	{
		for (int i = 0; i < this->atom_numbers; i = i + 1)
		{
			fscanf(fin, "%f %f %f",
				&this->velocity[i].x,
				&this->velocity[i].y,
				&this->velocity[i].z);
		}
	}
	if (irest == 0 || !has_vel)
	{
		for (int i = 0; i < this->atom_numbers; i = i + 1)
		{
			this->velocity[i].x = 0.0;
			this->velocity[i].y = 0.0;
			this->velocity[i].z = 0.0;
		}
	}
	fscanf(fin, "%f %f %f", &this->sys.box_length.x, &this->sys.box_length.y, &this->sys.box_length.z);
	controller.printf("        system size is %f %f %f\n", this->sys.box_length.x, this->sys.box_length.y, this->sys.box_length.z);
	cudaMemcpy(this->crd, this->coordinate, sizeof(VECTOR)*this->atom_numbers, cudaMemcpyHostToDevice);
	cudaMemcpy(this->vel, this->velocity, sizeof(VECTOR)*this->atom_numbers, cudaMemcpyHostToDevice);
	//in some bad rst7, the coordinates will be extremly bad, so need a full box map
	for (int i = 0; i < 10; i = i + 1)
	{
		Crd_Periodic_Map << <ceilf((float)this->atom_numbers / 32), 32 >> >
			(this->atom_numbers, this->crd, this->sys.box_length);
	}
	fclose(fin);
	controller.printf("    End reading AMBER rst7\n\n");
	
}

void MD_INFORMATION::trajectory_output::Append_Crd_Traj_File(FILE *fp)
{
	if (md_info->is_initialized)
	{
		md_info->Crd_Vel_Device_To_Host();
		if (fp == NULL)
		{
			fp = crd_traj;
		}
		fwrite(&md_info->coordinate[0].x, sizeof(VECTOR), md_info->atom_numbers, fp);
	}
}

void MD_INFORMATION::trajectory_output::Append_Box_Traj_File(FILE *fp)
{
	if (md_info->is_initialized)
	{
		if (fp == NULL)
		{
			fp = box_traj;
		}
		fprintf(fp, "%f %f %f %.0f %.0f %.0f\n", md_info->sys.box_length.x, md_info->sys.box_length.y, md_info->sys.box_length.z, 90.0f, 90.0f, 90.0f);
	}
}

void MD_INFORMATION::trajectory_output::Export_Restart_File(const char *rst7_name)
{
	if (!md_info->is_initialized)
		return;

	char filename[CHAR_LENGTH_MAX];
	if (rst7_name == NULL)
		strcpy(filename, restart_name);
	else
		strcpy(filename, rst7_name);
	md_info->Crd_Vel_Device_To_Host();
	if (amber_irest >= 0)
	{
		const char *sys_name = md_info->md_name;
		FILE *lin = NULL;
		Open_File_Safely(&lin, filename, "w");
		fprintf(lin, "%s\n", sys_name);
		fprintf(lin, "%8d %.3lf\n", md_info->atom_numbers, md_info->sys.Get_Current_Time());
		int s = 0;
		for (int i = 0; i < md_info->atom_numbers; i = i + 1)
		{
			fprintf(lin, "%12.7f%12.7f%12.7f", md_info->coordinate[i].x, md_info->coordinate[i].y, md_info->coordinate[i].z);
			s = s + 1;
			if (s == 2)
			{
				s = 0;
				fprintf(lin, "\n");
			}
		}
		if (s == 1)
		{
			s = 0;
			fprintf(lin, "\n");
		}
		for (int i = 0; i < md_info->atom_numbers; i = i + 1)
		{
			fprintf(lin, "%12.7f%12.7f%12.7f", md_info->velocity[i].x, md_info->velocity[i].y, md_info->velocity[i].z);
			s = s + 1;
			if (s == 2)
			{
				s = 0;
				fprintf(lin, "\n");
			}
		}
		if (s == 1)
		{
			s = 0;
			fprintf(lin, "\n");
		}
		fprintf(lin, "%12.7f%12.7f%12.7f", (float)md_info->sys.box_length.x, (float)md_info->sys.box_length.y, (float)md_info->sys.box_length.z);
		fprintf(lin, "%12.7f%12.7f%12.7f", (float)90., (float)90., (float)90.);
		fclose(lin);
	}
	else
	{
		FILE *lin = NULL;
		FILE *lin2 = NULL;
		char buffer[CHAR_LENGTH_MAX];
		sprintf(buffer, "%s_%s.txt", filename, "coordinate");
		Open_File_Safely(&lin, buffer, "w");
		sprintf(buffer, "%s_%s.txt", filename, "velocity");
		Open_File_Safely(&lin2, buffer, "w");
		fprintf(lin, "%d %.3lf\n", md_info->atom_numbers, md_info->sys.Get_Current_Time());
		fprintf(lin2, "%d %.3lf\n", md_info->atom_numbers, md_info->sys.Get_Current_Time());
		for (int i = 0; i < md_info->atom_numbers; i++)
		{
			fprintf(lin, "%12.7f %12.7f %12.7f\n", md_info->coordinate[i].x, md_info->coordinate[i].y, md_info->coordinate[i].z);
			fprintf(lin2, "%12.7f %12.7f %12.7f\n", md_info->velocity[i].x, md_info->velocity[i].y, md_info->velocity[i].z);
		}
		fprintf(lin, "%12.7f %12.7f %12.7f %12.7f %12.7f %12.7f", md_info->sys.box_length.x, md_info->sys.box_length.y, md_info->sys.box_length.z, 90.0f, 90.0f, 90.0f);
		fclose(lin);
		fclose(lin2);
	}
}


void MD_INFORMATION::Update_Volume(double factor)
{
	double f_inv = 1.0 / factor;

	sys.box_length = factor * sys.box_length;
	pbc.crd_to_uint_crd_cof = CONSTANT_UINT_MAX_FLOAT / sys.box_length;
	pbc.quarter_crd_to_uint_crd_cof = 0.25 * pbc.crd_to_uint_crd_cof;
	pbc.uint_dr_to_dr_cof = 1.0f / pbc.crd_to_uint_crd_cof;
	MD_Information_Crd_To_Uint_Crd();
}


float MD_INFORMATION::system_information::Get_Density()
{
	density = total_mass * 1e24f / 6.023e23f / Get_Volume();
	return density;
}

double MD_INFORMATION::system_information::Get_Current_Time()
{
	current_time = start_time + (double)dt_in_ps * steps;
	return current_time;
}

float MD_INFORMATION::system_information::Get_Volume()
{
	volume = box_length.x * box_length.y * box_length.z;
	return volume;
}

void MD_INFORMATION::MD_Information_Crd_To_Uint_Crd()
{
	Crd_To_Uint_Crd << <ceilf((float)this->atom_numbers / 128), 128 >> >
		(this->atom_numbers, pbc.quarter_crd_to_uint_crd_cof, crd, uint_crd);
}


void MD_INFORMATION::NVE_iteration::Leap_Frog()
{
	if (max_velocity <= 0)
	{
		MD_Iteration_Leap_Frog << <ceilf((float)md_info->atom_numbers / 128), 128 >> >
			(md_info->atom_numbers, md_info->vel, md_info->crd, md_info->frc, md_info->acc, md_info->d_mass_inverse, md_info->dt);
	}
	else
	{
		MD_Iteration_Leap_Frog_With_Max_Velocity << <ceilf((float)md_info->atom_numbers / 128), 128 >> >
			(md_info->atom_numbers, md_info->vel, md_info->crd, md_info->frc, md_info->acc, md_info->d_mass_inverse, md_info->dt, max_velocity);
	}
}

void MD_INFORMATION::MD_Information_Gradient_Descent()
{
	MD_Iteration_Gradient_Descent << <ceilf((float)this->atom_numbers / 128), 128 >> >
		(this->atom_numbers, this->crd, this->frc, dt * dt);
	Reset_List << <ceilf((float)this->atom_numbers * 3 / 128), 128 >> > (atom_numbers * 3, (float*)vel, 0);
}

void MD_INFORMATION::NVE_iteration::Velocity_Verlet_1()
{
	MD_Iteration_Speed_Verlet_1 << <ceilf((float)md_info->atom_numbers / 32), 32 >> >
		(md_info->atom_numbers, 0.5*md_info->dt, md_info->dt, md_info->acc, md_info->vel, md_info->crd, md_info->frc);
}

void MD_INFORMATION::NVE_iteration::Velocity_Verlet_2()
{
	if (max_velocity <= 0)
	{
		MD_Iteration_Speed_Verlet_2 << <ceilf((float)md_info->atom_numbers / 32), 32 >> >
			(md_info->atom_numbers, 0.5*md_info->dt, md_info->d_mass_inverse, md_info->frc, md_info->vel, md_info->acc);
	}
	else
	{
		MD_Iteration_Speed_Verlet_2_With_Max_Velocity << <ceilf((float)md_info->atom_numbers / 32), 32 >> >
			(md_info->atom_numbers, 0.5*md_info->dt, md_info->d_mass_inverse, md_info->frc, md_info->vel, md_info->acc, max_velocity);
	}
}


float MD_INFORMATION::system_information::Get_Total_Atom_Ek(int is_download)
{
	MD_Atom_Ek << <ceilf((float)md_info->atom_numbers / 32.), 32 >> >(md_info->atom_numbers, md_info->d_atom_ek, md_info->vel, md_info->d_mass);
	Sum_Of_List << <1, 1024 >> >(md_info->atom_numbers, md_info->d_atom_ek, d_sum_of_atom_ek);
	if (is_download)
	{
		cudaMemcpy(&h_sum_of_atom_ek, d_sum_of_atom_ek, sizeof(float), cudaMemcpyDeviceToHost);
		return h_sum_of_atom_ek;
	}
	else
	{
		return 0;
	}	
}

float MD_INFORMATION::system_information::Get_Atom_Temperature()
{
	h_temperature = Get_Total_Atom_Ek() * 2. / CONSTANT_kB / freedom;
	return h_temperature;
}

float MD_INFORMATION::residue_information::Get_Total_Residue_Ek(int is_download)
{
	MD_Residue_Ek << <ceilf((float)residue_numbers / 32.), 32 >> >
		(residue_numbers, d_res_start, d_res_end, res_ek_energy,
		md_info->vel, md_info->d_mass);
	Sum_Of_List << <1, 1024 >> >(residue_numbers, res_ek_energy, sigma_of_res_ek);
	if (is_download)
	{
		cudaMemcpy(h_sigma_of_res_ek, sigma_of_res_ek, sizeof(float), cudaMemcpyDeviceToHost);
		return h_sigma_of_res_ek[0];
	}
	else
	{
		return 0;
	}
}

float MD_INFORMATION::residue_information::Get_Residue_Temperature()
{
	h_temperature = Get_Total_Residue_Ek() * 2. / CONSTANT_kB / residue_numbers / 3;
	return h_temperature;
}

void MD_INFORMATION::residue_information::Residue_Crd_Map(VECTOR *no_wrap_crd, float scaler)
{
	Get_Center_Of_Mass << <20, 32 >> >(residue_numbers, d_res_start, d_res_end, no_wrap_crd, md_info->d_mass, d_mass_inverse, d_center_of_mass);
	Map_Center_Of_Mass << <20, { 32, 4 } >> >(residue_numbers, d_res_start, d_res_end, scaler, d_center_of_mass, md_info->sys.box_length, no_wrap_crd, md_info->crd);
}


void MD_INFORMATION::MD_Reset_Atom_Energy_And_Virial()
{
	need_potential = 0;
	Reset_List << <ceilf((float)atom_numbers / 1024.0f), 1024 >> >(atom_numbers, d_atom_energy, 0.0f);
	Reset_List << <1, 1 >> >(1, sys.d_potential, 0.0f);

	need_pressure = 0;
	Reset_List << <ceilf((float)atom_numbers / 1024.0f), 1024 >> >(atom_numbers, d_atom_virial, 0.0f);
	Reset_List << <1, 1 >> >(1, sys.d_virial, 0.0f);

	Reset_List << <ceilf((float)3 * atom_numbers / 1024.0f), 1024 >> >(3 * atom_numbers, (float*)frc, 0.0f);
}


void MD_INFORMATION::Calculate_Pressure_And_Potential_If_Needed(int is_download)
{
	if (need_pressure > 0)
	{
		sys.Get_Pressure(is_download);
	}
	if (need_potential > 0)
	{
		sys.Get_Potential(is_download);
	}
}

float MD_INFORMATION::system_information::Get_Pressure(int is_download)
{
	//计算动能
	MD_Atom_Ek << <ceilf((float)md_info->atom_numbers / 32.), 32 >> >(md_info->atom_numbers, md_info->d_atom_ek, md_info->vel, md_info->d_mass);
	Sum_Of_List << <1, 1024 >> >(md_info->atom_numbers, md_info->d_atom_ek, d_sum_of_atom_ek);

	//计算维里
	Add_Sum_List << <1, 1024 >> >(md_info->atom_numbers, md_info->d_atom_virial, d_virial);

	//合并起来
	Calculate_Pressure_Cuda << <1, 1 >> >(1.0 / Get_Volume(), d_sum_of_atom_ek, d_virial, d_pressure);

	if (is_download)
	{
		cudaMemcpy(&h_pressure, d_pressure, sizeof(float), cudaMemcpyDeviceToHost);
		return h_pressure;
	}
	else
	{
		return 0;
	}
}

float MD_INFORMATION::system_information::Get_Potential(int is_download)
{

	Add_Sum_List << <1, 1024 >> >(md_info->atom_numbers, md_info->d_atom_energy, d_potential);

	if (is_download)
	{
		cudaMemcpy(&h_potential, d_potential, sizeof(float), cudaMemcpyDeviceToHost);
		return h_potential;
	}
	else
	{
		return 0;
	}
}

void MD_INFORMATION::MD_Information_Frc_Device_To_Host()
{
	cudaMemcpy(this->force, this->frc, sizeof(VECTOR)*this->atom_numbers, cudaMemcpyDeviceToHost);
}

void MD_INFORMATION::MD_Information_Frc_Host_To_Device()
{
	cudaMemcpy(this->frc, this->force, sizeof(VECTOR)*this->atom_numbers, cudaMemcpyHostToDevice);
}

void MD_INFORMATION::Crd_Vel_Device_To_Host(int Do_Translation, int forced)
{
	if (output.current_crd_synchronized_step != sys.steps || forced)
	{	
		output.current_crd_synchronized_step = sys.steps;
		if (Do_Translation)
		{	
			cudaMemcpy(coordinate, crd, sizeof(VECTOR) * atom_numbers, cudaMemcpyDeviceToHost);
		}
		else
		{
			cudaMemcpy(this->coordinate, this->crd, sizeof(VECTOR)*this->atom_numbers, cudaMemcpyDeviceToHost);
		}
		cudaMemcpy(this->velocity, this->vel, sizeof(VECTOR)*this->atom_numbers, cudaMemcpyDeviceToHost);
	}
}

void MD_INFORMATION::Clear()
{

}
