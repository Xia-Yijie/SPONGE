#include "TI_core.cuh"

#define TRAJ_COMMAND "crd"
#define TRAJ_DEFAULT_FILENAME "mdcrd.dat"
//#define VEL_COMMAND "vel"
//#define VEL_DEFAULT_FILENAME "mdvel.dat"
#define BOX_COMMAND "box"
#define BOX_DEFAULT_FILENAME "box.txt"


void TI_CORE::non_bond_information::Initial(CONTROLLER *controller, TI_CORE *TI_core )
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
		int toscan = fscanf(fp, "%d %d", &atom_numbers, &excluded_atom_numbers);
		if (TI_core->atom_numbers > 0 && TI_core->atom_numbers != atom_numbers)
		{
			controller->printf("        Error: atom_numbers is not equal: %d %d\n", TI_core->atom_numbers, atom_numbers);
			getchar();
			exit(1);
		}
		else if (TI_core->atom_numbers == 0)
		{
			TI_core->atom_numbers = atom_numbers;
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
			toscan = fscanf(fp, "%d", &h_excluded_numbers[i]);
			h_excluded_list_start[i] = count;
			for (int j = 0; j < h_excluded_numbers[i]; j++)
			{
				toscan = fscanf(fp, "%d", &h_excluded_list[count]);
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
				char *toget = fgets(temps, CHAR_LENGTH_MAX, parm);

				int atom_numbers = 0;
				int toscan = fscanf(parm, "%d\n", &atom_numbers);
				if (TI_core->atom_numbers > 0 && TI_core->atom_numbers != atom_numbers)
				{
					controller->printf("        Error: atom_numbers is not equal: %d %d\n", TI_core->atom_numbers, atom_numbers);
					getchar();
					exit(1);
				}
				else if (TI_core->atom_numbers == 0)
				{
					TI_core->atom_numbers = atom_numbers;
				}
				Cuda_Malloc_Safely((void**)&d_excluded_list_start, sizeof(int)*atom_numbers);
				Cuda_Malloc_Safely((void**)&d_excluded_numbers, sizeof(int)*atom_numbers);

				Malloc_Safely((void**)&h_excluded_list_start, sizeof(int)*atom_numbers);
				Malloc_Safely((void**)&h_excluded_numbers, sizeof(int)*atom_numbers);
				for (int i = 0; i < 9; i = i + 1)
				{
					toscan = fscanf(parm, "%d\n", &excluded_atom_numbers);
				}
				toscan = fscanf(parm, "%d\n", &excluded_atom_numbers);
				controller->printf("        excluded list total length is %d\n", excluded_atom_numbers);

				Cuda_Malloc_Safely((void**)&d_excluded_list, sizeof(int)*excluded_atom_numbers);
				Malloc_Safely((void**)&h_excluded_list, sizeof(int)*excluded_atom_numbers);
			}

			//read atom_excluded_number for every atom
			if (strcmp(temp_first_str, "%FLAG") == 0
				&& strcmp(temp_second_str, "NUMBER_EXCLUDED_ATOMS") == 0)
			{
				char *toget = fgets(temps, CHAR_LENGTH_MAX, parm);
				for (int i = 0; i<TI_core->atom_numbers; i = i + 1)
				{
					int toscan = fscanf(parm, "%d\n", &h_excluded_numbers[i]);
				}
			}
			//read every atom's excluded atom list
			if (strcmp(temp_first_str, "%FLAG") == 0
				&& strcmp(temp_second_str, "EXCLUDED_ATOMS_LIST") == 0)
			{
				int count = 0;
				//int none_count = 0;
				int lin = 0;
				char *toget = fgets(temps, CHAR_LENGTH_MAX, parm);
				for (int i = 0; i<TI_core->atom_numbers; i = i + 1)
				{
					h_excluded_list_start[i] = count;
					for (int j = 0; j<h_excluded_numbers[i]; j = j + 1)
					{
						int toscan = fscanf(parm, "%d\n", &lin);
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

		cudaMemcpy(d_excluded_list_start, h_excluded_list_start, sizeof(int)*TI_core->atom_numbers, cudaMemcpyHostToDevice);
		cudaMemcpy(d_excluded_numbers, h_excluded_numbers, sizeof(int)*TI_core->atom_numbers, cudaMemcpyHostToDevice);
		cudaMemcpy(d_excluded_list, h_excluded_list, sizeof(int)*excluded_atom_numbers, cudaMemcpyHostToDevice);
		controller->printf("    End reading excluded list from AMBER file\n\n");
		fclose(parm);
	}
	else
	{
		int atom_numbers = TI_core->atom_numbers;
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

void TI_CORE::periodic_box_condition_information::Update_Volume(VECTOR box_length)
{
	crd_to_uint_crd_cof = CONSTANT_UINT_MAX_FLOAT / box_length;
	quarter_crd_to_uint_crd_cof = 0.25 * crd_to_uint_crd_cof;
	uint_dr_to_dr_cof = 1.0f / crd_to_uint_crd_cof;
}

void TI_CORE::trajectory_input::Initial(CONTROLLER *controller, TI_CORE * TI_core)
{
	this->TI_core = TI_core;
	if (controller[0].Command_Exist("frame_numbers"))
	{
		frame_numbers = atoi(controller[0].Command("frame_numbers"));
	}
	else
	{
		printf("	warning: missing value of frame numbers, set to default 1000.\n");
		frame_numbers = 1000;
	}
	current_frame = 0;
	bytes_per_frame = TI_core->atom_numbers * 3 * sizeof(float);
	if (controller[0].Command_Exist(TRAJ_COMMAND))
	{
		Open_File_Safely(&crd_traj, controller[0].Command(TRAJ_COMMAND), "rb");
	}
	else
	{
		printf("	Error: missing trajectory file.\n");
		getchar();
		exit(1);
	}
	if (controller[0].Command_Exist(BOX_COMMAND))
	{
		Open_File_Safely(&box_traj, controller[0].Command(BOX_COMMAND), "r");
	}
	else
	{
		printf("	Error: missing box trajectory file.\n");
		getchar();
		exit(1);
	}

	/*if (controller[0].Command_Exist("mass_pertubated"))
	{
		TI_core->mass_pertubated = atoi(controller[0].Command("mass_pertubated"));
	}
	else
	{
		printf("	missing value of mass_pertubated, set to defualt 0.");
		TI_core->mass_pertubated = 0;
	}
	if (TI_core->mass_pertubated != 0)
	{
		if (controller[0].Command_Exist(VEL_COMMAND))
		{
			Open_File_Safely(&vel_traj, controller[0].Command(VEL_COMMAND), "rb");
		}
		else
		{
			printf("	Error: missing velocity trajectory file.\n");
			getchar();
			exit(1);
		}
	}*/
}

void TI_CORE::Initial(CONTROLLER *controller)
{
	controller[0].printf("START INITIALZING TI CORE:\n");
	
	if (controller[0].Command_Exist("atom_numbers"))
	{
		atom_numbers = atoi(controller[0].Command("atom_numbers"));
	}
	else
	{
		printf("	Error: missing value of atom numbers.\n");
		getchar();
		exit(1);
	}

	box_length.x = box_length.y = box_length.z = 1.0;
	last_box_length.x = last_box_length.y = last_box_length.z = 1.0;
	volume_change_factor = 0.0;
	box_angle.x = box_angle.y = box_angle.z = 0.0;
	if (controller[0].Command_Exist("charge_pertubated"))
	{
		charge_pertubated = atoi(controller[0].Command("charge_pertubated"));
	}
	else
	{
		printf("	Warning: missing value of charge pertubated, set to default 0.\n");
		charge_pertubated = 0;
	}
	if (controller[0].Command_Exist("bond_pertubated"))
	{
		bond_pertubated = atoi(controller[0].Command("bond_pertubated"));
	}
	else
	{
		printf("	Warning: missing value of bond pertubated, set to default 0.\n");
		bond_pertubated = 0;
	}
	if (controller[0].Command_Exist("angle_pertubated"))
	{
		angle_pertubated = atoi(controller[0].Command("angle_pertubated"));
	}
	else
	{
		printf("	Warning: missing value of angle pertubated, set to default 0.\n");
		angle_pertubated = 0;
	}
	if (controller[0].Command_Exist("dihedral_pertubated"))
	{
		dihedral_pertubated = atoi(controller[0].Command("dihedral_pertuabted"));
	}
	else
	{
		printf("	Warning: missing value of dihedral pertubated, set to default 0.\n");
		dihedral_pertubated = 0;
	}
	if (controller[0].Command_Exist("nb14_pertubated"))
	{
		nb14_pertubated = atoi(controller[0].Command("nb14_pertubated"));
	}
	else
	{
		printf("	Warning: missing value of nb14 pertubated, set to default 0.\n");
		nb14_pertubated = 0;
	}
	if (controller[0].Command_Exist("lj_pertubated"))
	{
		lj_pertubated = atoi(controller[0].Command("lj_pertubated"));
	}
	else
	{
		printf("	Warning: missing value of lj pertubated, set to default 0.\n");
		lj_pertubated = 0;
	}

	Malloc_Safely((void**)&h_charge, sizeof(float) * atom_numbers);
	Malloc_Safely((void**)&h_charge_A, sizeof(float) * atom_numbers);
	Malloc_Safely((void**)&h_charge_B, sizeof(float) * atom_numbers);
	Malloc_Safely((void**)&h_charge_B_A, sizeof(float) * atom_numbers);
	Malloc_Safely((void**)&h_subsys_division, sizeof(int) * atom_numbers);
	Cuda_Malloc_Safely((void**)&d_charge, sizeof(float) * atom_numbers);
	Cuda_Malloc_Safely((void**)&d_charge_B_A, sizeof(float) * atom_numbers);
	Cuda_Malloc_Safely((void**)&d_subsys_division, sizeof(int) * atom_numbers);
	
	need_nbl = ((charge_pertubated > 0) || (lj_pertubated) > 0);
	
	if (charge_pertubated > 0)
	{
		if (controller[0].Command_Exist("chargeA_in_file") && controller[0].Command_Exist("chargeB_in_file"))
		{
			controller[0].printf("    Start reading chargeA:\n");
			int atom_numbers_in_file = 0;
			FILE * fp = NULL;
			Open_File_Safely(&fp, controller[0].Command("chargeA_in_file"), "r");
			char lin[CHAR_LENGTH_MAX];
			char *toget = fgets(lin, CHAR_LENGTH_MAX, fp);
			int scanf_ret = sscanf(lin, "%d", &atom_numbers_in_file);
			if (atom_numbers != atom_numbers_in_file)
			{
				controller->printf("        Error: atom_numbers is not equal: %d %d\n", atom_numbers, atom_numbers_in_file);
				getchar();
				exit(1);
			}
			for (int i = 0; i < atom_numbers; ++i)
			{
				scanf_ret = fscanf(fp, "%f", &h_charge_A[i]);
			}
			fclose(fp);
			controller[0].printf("    End reading chargeA\n\n");

			controller[0].printf("    Start reading chargeB:\n");
			Open_File_Safely(&fp, controller[0].Command("chargeB_in_file"), "r");
			toget = fgets(lin, CHAR_LENGTH_MAX, fp);
			scanf_ret = sscanf(lin, "%d", &atom_numbers_in_file);
			if (atom_numbers != atom_numbers_in_file)
			{
				controller->printf("        Error: atom_numbers is not equal: %d %d\n", atom_numbers, atom_numbers_in_file);
				getchar();
				exit(1);
			}
			for (int i = 0; i < atom_numbers; ++i)
			{
				scanf_ret = fscanf(fp, "%f", &h_charge_B[i]);
			}
			fclose(fp);
			controller[0].printf("    End reading chargeB\n\n");
			for (int i = 0; i < atom_numbers; ++i)
			{
				h_charge_B_A[i] = h_charge_B[i] - h_charge_A[i];
			}
			cudaMemcpy(d_charge_B_A, h_charge_B_A, sizeof(float)*atom_numbers, cudaMemcpyHostToDevice);
		}
		else
		{
			printf("	Error: missing value of charge A and charge B, These value must be given in TI mode if charge is pertubated.\n");
			getchar();
			exit(1);
		}
	}

	if (controller[0].Command_Exist("charge_in_file"))
	{
		FILE *fp = NULL;
		controller->printf("    Start reading charge:\n");
		Open_File_Safely(&fp, controller[0].Command("charge_in_file"), "r");
		int atom_numbers = 0;
		char lin[CHAR_LENGTH_MAX];
		char *toget = fgets(lin, CHAR_LENGTH_MAX, fp);
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
		for (int i = 0; i < atom_numbers; i++)
		{
			scanf_ret = fscanf(fp, "%f", &h_charge[i]);
		}
		controller->printf("    End reading charge\n\n");
		fclose(fp);
	}
	else if (atom_numbers > 0)
	{
		controller[0].printf("    charge is set to 0 as default\n");
		for (int i = 0; i < atom_numbers; i++)
		{
			h_charge[i] = 0;
		}
	}
	cudaMemcpy(d_charge, h_charge, sizeof(float)*atom_numbers, cudaMemcpyHostToDevice);


	if (controller[0].Command_Exist("subsys_division_in_file"))
	{
		FILE *fp = NULL;
		controller->printf("    Start reading subsystem division information:\n");
		Open_File_Safely(&fp, controller[0].Command("subsys_division_in_file"), "r");
		int atom_numbers = 0;
		char lin[CHAR_LENGTH_MAX];
		char *toget = fgets(lin, CHAR_LENGTH_MAX, fp);
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
		for (int i = 0; i < atom_numbers; i++)
		{
			scanf_ret = fscanf(fp, "%d", &h_subsys_division[i]);
		}
		controller->printf("    End reading subsystem information\n\n");
		fclose(fp);
	}
	else if (atom_numbers > 0)
	{
		controller[0].printf("    subsystem mask is set to 0 as default\n");
		for (int i = 0; i < atom_numbers; i++)
		{
			h_subsys_division[i] = 0;
		}
	}
	cudaMemcpy(d_subsys_division, h_subsys_division, sizeof(float)*atom_numbers, cudaMemcpyHostToDevice);

	//Malloc_Safely((void**)&velocity, sizeof(VECTOR) * atom_numbers);
	Malloc_Safely((void**)&coordinate, sizeof(VECTOR) * atom_numbers);
	//Cuda_Malloc_Safely((void**)&vel, sizeof(VECTOR) * atom_numbers);
	Cuda_Malloc_Safely((void**)&crd, sizeof(VECTOR) * atom_numbers);
	Cuda_Malloc_Safely((void**)&uint_crd, sizeof(UNSIGNED_INT_VECTOR)*atom_numbers);

	nb.Initial(controller, this);
	input.Initial(controller, this);

	controller[0].Step_Print_Initial("frame", "%d");
	controller[0].Step_Print_Initial("dH_dlambda", "%.2f");
	Read_Next_Frame();

	printf("END INITIALZING TI CORE\n\n");
}

void TI_CORE::TI_Core_Crd_To_Uint_Crd()
{
	Crd_To_Uint_Crd << <ceilf((float)this->atom_numbers / 128), 128 >> >
		(this->atom_numbers, pbc.quarter_crd_to_uint_crd_cof, crd, uint_crd);
}

void TI_CORE::Read_Next_Frame()
{
	size_t toread = fread(coordinate, sizeof(VECTOR), atom_numbers, input.crd_traj);
	cudaMemcpy(crd, coordinate, sizeof(VECTOR) *atom_numbers, cudaMemcpyHostToDevice);
	last_box_length.x = box_length.x;
	last_box_length.y = box_length.y;
	last_box_length.z = box_length.z;
	int toscan = fscanf(input.box_traj, "%f %f %f %f %f %f", &box_length.x, &box_length.y, &box_length.z, &box_angle.x, &box_angle.y, &box_angle.z);
	volume_change_factor = box_length.x / last_box_length.x;
	/*if (mass_pertubated != 0)
	{
		fread(velocity, sizeof(VECTOR), atom_numbers, input.vel_traj);
		cudaMemcpy(vel, velocity, sizeof(VECTOR) * atom_numbers, cudaMemcpyHostToDevice);
	}*/
	pbc.Update_Volume(box_length);
	TI_Core_Crd_To_Uint_Crd();
}


void TI_CORE::Clear()
{
	//free(velocity);
	free(coordinate);
	//cudaFree(vel);
	cudaFree(crd);
	cudaFree(uint_crd);

	free(h_charge_A);
	free(h_charge_B);
	free(h_charge);
	free(h_charge_B_A);
	cudaFree(d_charge);
	cudaFree(d_charge_B_A);
	free(h_subsys_division);
	cudaFree(d_subsys_division);
	
	//velocity = NULL;
	coordinate = NULL;
	//vel = NULL;
	crd = NULL;
	uint_crd = NULL;
	h_charge_A = NULL;
	h_charge_B = NULL;
	h_charge = NULL;
	h_charge_B_A = NULL;
	d_charge = NULL;
	d_charge_B_A = NULL;
	h_subsys_division = NULL;
	d_subsys_division = NULL;

	free(nb.h_excluded_list);
	free(nb.h_excluded_numbers);
	free(nb.h_excluded_list_start);
	cudaFree(nb.d_excluded_list);
	cudaFree(nb.d_excluded_numbers);
	cudaFree(nb.d_excluded_list_start);
	nb.h_excluded_list = NULL;
	nb.h_excluded_numbers = NULL;
	nb.h_excluded_list_start = NULL;
	nb.d_excluded_list_start = NULL;
	nb.d_excluded_numbers = NULL;
	nb.d_excluded_list = NULL;

	fclose(input.crd_traj);
	/*if (mass_pertubated)
		fclose(input.vel_traj);*/
	fclose(input.box_traj);
}

void TI_CORE::dH_dlambda_data::Sum_One_Frame()
{
	dH_dlambda_current_frame = bond_soft_dH_dlambda + lj_soft_dH_dlambda + coul_direct_dH_dlambda +lj_soft_long_range_correction + pme_dH_dlambda /*+ kinetic_dH_dlambda*/ + (bondB_ene - bondA_ene) + (angleB_ene - angleA_ene) + (dihedralB_ene - dihedralA_ene) + (nb14B_EE_ene - nb14A_EE_ene) + (nb14B_LJ_ene - nb14A_LJ_ene);
	
	total_dH_dlambda += dH_dlambda_current_frame;
}

void TI_CORE::TI_Core_Crd_Device_To_Host()
{
	cudaMemcpy(coordinate, crd, sizeof(VECTOR) * atom_numbers, cudaMemcpyDeviceToHost);
}

void TI_CORE::Print_dH_dlambda_Average_To_Screen_And_Result_File(FILE * fcresult)
{
	data.average_dH_dlambda = data.total_dH_dlambda/input.frame_numbers;
	fprintf(stdout, "Ensemble Average <dH/dlambda>: %.6f\n", data.average_dH_dlambda);
	fprintf(fcresult, "%.6f\n", data.average_dH_dlambda);
}


