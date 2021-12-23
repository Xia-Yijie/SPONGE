#include "constrain.cuh"


void CONSTRAIN::Initial_Constrain(CONTROLLER *controller, const int atom_numbers, const float dt, const VECTOR box_length, const float exp_gamma, const int is_Minimization, float *atom_mass, int *system_freedom)
{
	//�Ӵ���Ĳ������ƻ�����Ϣ
	this->atom_numbers = atom_numbers;
	this->dt = dt;
	this->dt_inverse = 1.0 / dt;
	this->quarter_crd_to_uint_crd_cof = 0.25 * CONSTANT_UINT_MAX_FLOAT / box_length;
	this->uint_dr_to_dr_cof = 1.0f / CONSTANT_UINT_MAX_FLOAT * box_length;
	this->volume = box_length.x * box_length.y * box_length.z;


	//ȷ��ʹ�õ���֮����ԡ��������exp_gamma
	this->v_factor = exp_gamma;
	this->x_factor = 0.5*(1. + exp_gamma);

	if (is_Minimization)
	{
		this->v_factor = 0.0f;
	}

	int extra_numbers = 0;
	FILE *fp = NULL;
	//���ļ���һ����ȷ��constrain������Ϊ�����ڴ���׼��
	if (controller[0].Command_Exist(this->module_name, "in_file"))
	{
		Open_File_Safely(&fp, controller[0].Command(this->module_name, "in_file"), "r");
		int scanf_ret = fscanf(fp, "%d", &extra_numbers);
	}

	constrain_pair_numbers = bond_constrain_pair_numbers + angle_constrain_pair_numbers + extra_numbers;
	system_freedom[0] -= constrain_pair_numbers;
	controller[0].printf("    constrain pair number is %d\n", constrain_pair_numbers);

	Malloc_Safely((void**)&h_constrain_pair, sizeof(CONSTRAIN_PAIR)*constrain_pair_numbers);
	Cuda_Malloc_Safely((void**)&constrain_pair, sizeof(CONSTRAIN_PAIR)*constrain_pair_numbers);
	for (int i = 0; i < bond_constrain_pair_numbers; i = i + 1)
	{
		h_constrain_pair[i] = h_bond_pair[i];
		h_constrain_pair[i].constrain_k = h_constrain_pair[i].constrain_k / this->x_factor;
	}
	for (int i = 0; i < angle_constrain_pair_numbers; i = i + 1)
	{
		h_constrain_pair[i + bond_constrain_pair_numbers] = h_angle_pair[i];
		h_constrain_pair[i + bond_constrain_pair_numbers].constrain_k = h_constrain_pair[i + bond_constrain_pair_numbers].constrain_k / this->x_factor;
	}
	//���ļ�����
	if (fp != NULL)
	{
		int atom_i, atom_j;
		int count = bond_constrain_pair_numbers + angle_constrain_pair_numbers;
		for (int i = 0; i < extra_numbers; i = i + 1)
		{
			int scanf_ret = fscanf(fp, "%d %d %f", &atom_i, &atom_j, &h_constrain_pair[count].constant_r);
			h_constrain_pair[count].atom_i_serial = atom_i;
			h_constrain_pair[count].atom_j_serial = atom_j;
			h_constrain_pair[count].constrain_k = atom_mass[atom_i] * atom_mass[atom_j] / (atom_mass[atom_i] + atom_mass[atom_j]) / this->x_factor;
			count += 1;
		}
		fclose(fp);
		fp = NULL;
	}

	//�ϴ�GPU
	cudaMemcpy(constrain_pair, h_constrain_pair, sizeof(CONSTRAIN_PAIR)*constrain_pair_numbers, cudaMemcpyHostToDevice);


	//��ճ�ʼ��ʱʹ�õ���ʱ����
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

	if (is_initialized && !is_controller_printf_initialized)
	{
		is_controller_printf_initialized = 1;
		controller[0].printf("    structure last modify date is %d\n", last_modify_date);
	}
	controller[0].printf("END INITIALIZING CONSTRAIN\n\n");
	is_initialized = 1;
}


void CONSTRAIN::Add_HBond_To_Constrain_Pair
(CONTROLLER *controller, const int bond_numbers, const int *atom_a, const int *atom_b, const float *bond_r,
const float *atom_mass, char *module_name)
{
	controller[0].printf("START INITIALIZING CONSTRAIN:\n");
	if (module_name == NULL)
	{
		strcpy(this->module_name, "constrain");
	}
	else
	{
		strcpy(this->module_name, module_name);
	}
	constrain_mass = 3.3f;
	if (controller[0].Command_Exist(this->module_name,"in_file"))
		constrain_mass = 0.0f;
	if (controller[0].Command_Exist(this->module_name, "mass"))
		constrain_mass = atof(controller[0].Command(this->module_name, "mass"));
	//Ԥ�ȷ���һ���㹻���CONSTRAIN_PAIR������ʱ�洢
	Malloc_Safely((void**)&h_bond_pair, sizeof(CONSTRAIN_PAIR)*bond_numbers);
	int s = 0;
	float mass_a, mass_b;
	for (int i = 0; i < bond_numbers; i = i + 1)
	{
		mass_a = atom_mass[atom_a[i]];
		mass_b = atom_mass[atom_b[i]];
		if ((mass_a <constrain_mass && mass_a > 0) || (mass_b <constrain_mass && mass_b > 0))//����Hԭ�ӵ�bond
		{
			h_bond_pair[s].atom_i_serial = atom_a[i];
			h_bond_pair[s].atom_j_serial = atom_b[i];
			h_bond_pair[s].constant_r = bond_r[i];
			h_bond_pair[s].constrain_k = atom_mass[atom_a[i]] * atom_mass[atom_b[i]] / (atom_mass[atom_a[i]] + atom_mass[atom_b[i]]);
			s = s + 1;
		}
	}
	bond_constrain_pair_numbers = s;

	//����ʹ���߲����ڵ���Add_HBond_To_Constrain_Pair��Add_HAngle_To_Constrain_Pair��;�ͷ�atom_a��ָ��ָ��Ŀռ�
	bond_info.bond_numbers = bond_numbers;
	bond_info.atom_a = atom_a;
	bond_info.atom_b = atom_b;
	bond_info.bond_r = bond_r;
}
void CONSTRAIN::Add_HAngle_To_Constrain_Pair
(CONTROLLER *controller, const int angle_numbers, const int *atom_a, const int *atom_b, const int *atom_c,
const float *angle_theta, const float *atom_mass)
{
	int temp = 0;
	if (controller->Command_Exist(this->module_name, "angle") && atoi(controller->Command(this->module_name, "angle")))
	{
		temp = angle_numbers;
	}
	
	//Ĭ����Ϊ�Ѿ�������Add_HBond_To_Constrain_Pair

	//Ԥ�ȷ���һ���㹻���CONSTRAIN_PAIR������ʱ�洢
	Malloc_Safely((void**)&h_angle_pair, sizeof(CONSTRAIN_PAIR)*angle_numbers*2);
	int s = 0;
	float mass_a, mass_c;
	for (int i = 0; i < temp; i = i + 1)
	{
		mass_a = atom_mass[atom_a[i]];
		mass_c = atom_mass[atom_c[i]];
		if ((mass_a <constrain_mass && mass_a > 0) || (mass_c <constrain_mass && mass_c > 0))//����Hԭ�ӵ�angle,������ԭ�Ӳ����ڽ����ġ�
		{
			h_angle_pair[s].atom_i_serial = atom_a[i];//�̶�angle���˵�������
			h_angle_pair[s].atom_j_serial = atom_c[i];

			float rab=0., rbc=0.;
			for (int j = 0; j < bond_info.bond_numbers; j = j + 1)
			{
				//�ҵ�a��bԭ�ӵ�ƽ�����
				if ((bond_info.atom_a[j] == atom_a[i] && bond_info.atom_b[j] == atom_b[i])
					|| (bond_info.atom_a[j] == atom_b[i] && bond_info.atom_b[j] == atom_a[i]))
				{
					rab = bond_info.bond_r[j];
				}

				//�ҵ�b��cԭ�ӵ�ƽ�����
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

			//�������Ҷ���õ�ƽ���ac��������constrain
			h_angle_pair[s].constant_r = sqrtf(rab*rab + rbc*rbc - 2.*rab*rbc*cosf(angle_theta[i]));
			h_angle_pair[s].constrain_k = atom_mass[atom_a[i]] * atom_mass[atom_c[i]] / (atom_mass[atom_a[i]] + atom_mass[atom_c[i]]);

			s = s + 1;

		}
	}
	angle_constrain_pair_numbers = s;
}

void CONSTRAIN::Update_Volume(VECTOR box_length)
{
	if (is_initialized)
	{
		quarter_crd_to_uint_crd_cof = 0.25 * CONSTANT_UINT_MAX_FLOAT / box_length;
		uint_dr_to_dr_cof = 1.0f / CONSTANT_UINT_MAX_FLOAT * box_length;
		volume = box_length.x * box_length.y * box_length.z;
	}
}

void CONSTRAIN::Clear()
{
	if (is_initialized)
	{
		is_initialized = 0;

		cudaFree(constrain_pair);
		constrain_pair = NULL;

		free(h_constrain_pair);
		h_constrain_pair = NULL;
	}
}