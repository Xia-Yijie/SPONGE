#include "MC_barostat.cuh"

void MC_BAROSTAT_INFORMATION::Volume_Change_Attempt(VECTOR boxlength)
{
	double nrand = ((double)2.0 * rand() / RAND_MAX - 1.0);
	//double nrand = -1.0;
	DeltaV = nrand * DeltaV_max;
	double V = boxlength.x * boxlength.y * boxlength.z;
	newV = V + DeltaV;
	VDevided = newV / V;
	crd_scale_factor = cbrt(VDevided);
	//final_term = p0 * DeltaV - N_Beta_Inverse * logf(VDevided);
}


int MC_BAROSTAT_INFORMATION::Check_MC_Barostat_Accept()
{
	total_count += 1;
	if ( (float) rand() / RAND_MAX < accept_possibility)
	{
		reject = 0;
		accep_count += 1;
	}
	else
	{
		reject = 1;
	}
	return reject;
}

void MC_BAROSTAT_INFORMATION::Initial(CONTROLLER *controller, int atom_numbers, 
	float target_pressure, VECTOR boxlength, int res_is_initialized, char *module_name)
{
	controller->printf("START INITIALIZING MC BAROSTAT:\n");
	if (module_name == NULL)
	{
		strcpy(this->module_name, "mc_baro");
	}
	else
	{
		strcpy(this->module_name, module_name);
	}
	controller->printf("    The target pressure is %.2f bar\n", target_pressure*CONSTANT_PRES_CONVERTION);
	V0 = boxlength.x * boxlength.y * boxlength.z;
	newV = V0;
	float mc_baro_initial_ratio = 0.01;
	if (controller[0].Command_Exist(this->module_name, "initial_ratio"))
		mc_baro_initial_ratio = atof(controller[0].Command(this->module_name, "initial_ratio"));
	DeltaV_max = mc_baro_initial_ratio * V0;
	controller->printf("    The initial max volume to change is %f A^3\n", DeltaV_max);
	
	update_interval = 100;
	if (controller[0].Command_Exist(this->module_name, "update_interval"))
		update_interval = atoi(controller[0].Command(this->module_name, "update_interval"));
	controller->printf("    The update_interval is %d\n", update_interval);

	check_interval = 20;
	if (controller[0].Command_Exist(this->module_name, "check_interval"))
		check_interval = atoi(controller[0].Command(this->module_name, "check_interval"));
	controller->printf("    The check_interval is %d\n", check_interval);

	scale_coordinate_by_residue = res_is_initialized;
	if (controller[0].Command_Exist(this->module_name, "residue_scale"))
		scale_coordinate_by_residue = atoi(controller[0].Command(this->module_name, "residue_scale"));
	if (scale_coordinate_by_residue == 1 && res_is_initialized == 0)
	{
		controller->printf("    Warning: The residue is not initialized, so can not use residue scale mode. Atom scale mode is set instead.\n", update_interval);
		scale_coordinate_by_residue = 0;
	}
	controller->printf("    The residue_scale is %d\n", scale_coordinate_by_residue);

	system_reinitializing_count = 0;

	accept_rate_low = 30.0;
	if (controller[0].Command_Exist(this->module_name, "accept_rate_low"))
		accept_rate_low = atoi(controller[0].Command(this->module_name, "accept_rate_low"));
	controller->printf("    The lowest accept rate is %.2f%%\n", accept_rate_low);

	accept_rate_high = 40.0;
	if (controller[0].Command_Exist(this->module_name, "accept_rate_high"))
		accept_rate_high = atoi(controller[0].Command(this->module_name, "accept_rate_high"));
	controller->printf("    The highest accept rate is %.2f%%\n", accept_rate_high);


	Cuda_Malloc_Safely((void**)&frc_backup, sizeof(VECTOR)*atom_numbers);
	Cuda_Malloc_Safely((void**)&crd_backup, sizeof(VECTOR)*atom_numbers);
	is_initialized = 1;
	if (is_initialized && !is_controller_printf_initialized)
	{
		controller->Step_Print_Initial("density", "%.4f");
		is_controller_printf_initialized = 1;
		controller[0].printf("    structure last modify date is %d\n", last_modify_date);
	}

	controller[0].printf("END INITIALIZING MC BAROSTAT\n\n");
}

void MC_BAROSTAT_INFORMATION::Delta_V_Max_Update()
{
	if (total_count % check_interval == 0)
	{
		accept_rate = 100.0 * accep_count / total_count;

		if (accept_rate < accept_rate_low)
		{
			total_count = 0;
			accep_count = 0;
			DeltaV_max *= 0.9;
		}
		if (accept_rate > accept_rate_high)
		{
			total_count = 0;
			accep_count = 0;
			DeltaV_max *= 1.1;
		}
	}
}


void MC_BAROSTAT_INFORMATION::Ask_For_Calculate_Potential(int steps, int *need_potential)
{
	if (is_initialized && steps % update_interval == 0)
	{
		*need_potential = 1;
	}
}