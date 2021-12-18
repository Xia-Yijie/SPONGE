/*
* Copyright 2021 Gao's lab, Peking University, CCME. All rights reserved.
*
* NOTICE TO LICENSEE:
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
* http://www.apache.org/licenses/LICENSE-2.0
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/

#ifndef NB14_CUH
#define NB14_CUH
#include "../common.cuh"
#include "../control.cuh"

//用于计算LJ_Force时使用的坐标和记录的原子LJ种类序号与原子电荷
#ifndef UINT_VECTOR_LJ_TYPE_DEFINE
#define UINT_VECTOR_LJ_TYPE_DEFINE
struct UINT_VECTOR_LJ_TYPE
{
	unsigned int uint_x;
	unsigned int uint_y;
	unsigned int uint_z;
	int LJ_type;
	float charge;
};
__device__ __host__ VECTOR Get_Periodic_Displacement(const UINT_VECTOR_LJ_TYPE uvec_a, const UINT_VECTOR_LJ_TYPE uvec_b, const VECTOR scaler);
__global__ void Copy_LJ_Type_To_New_Crd(const int atom_numbers, UINT_VECTOR_LJ_TYPE *new_crd, const int *LJ_type);
__global__ void Copy_Crd_And_Charge_To_New_Crd(const int atom_numbers, const UNSIGNED_INT_VECTOR *crd, UINT_VECTOR_LJ_TYPE *new_crd, const float *charge);
__global__ void Copy_Crd_To_New_Crd(const int atom_numbers, const UNSIGNED_INT_VECTOR *crd, UINT_VECTOR_LJ_TYPE *new_crd);
#endif


struct NON_BOND_14
{
	char module_name[CHAR_LENGTH_MAX];
	int is_initialized = 0;
	int is_controller_printf_initialized = 0;
	int last_modify_date = 20210408;

	//r = ab原子的距离
	//E_lj_energy = lj_scale_factor * (lj_A/12 * r^-12 - lj_B/6 * r^-6) 
	//E_cf_energy = cf_scale_factor * charge_a * charge_b / r
	//lj_A、lj_B、charge从外部传入，lj_A、lj_B参考LJ，charge参考md_core
	int nb14_numbers = 0;
	int *h_atom_a = NULL;
	int *h_atom_b = NULL;
	int *d_atom_a = NULL;
	int *d_atom_b = NULL;
	float *h_lj_scale_factor = NULL;
	float *d_lj_scale_factor = NULL;
	float *h_cf_scale_factor = NULL;
	float *d_cf_scale_factor = NULL;

	float *d_nb14_energy = NULL;
	float *d_nb14_cf_energy_sum = NULL;
	float *d_nb14_lj_energy_sum = NULL;
	float h_nb14_cf_energy_sum = 0;
	float h_nb14_lj_energy_sum = 0;

	int threads_per_block = 128;

	void Initial(CONTROLLER *controller, char *module_name = NULL);
	void Clear();
	void Memory_Allocate();
	void Read_Information_From_AMBERFILE(const char *file_name, CONTROLLER controller);
	void Parameter_Host_To_Device();

	void Energy_Device_To_Host();

	/*-----------------------------------------------------------------------------------------
	下面的函数是普通md的需求
	------------------------------------------------------------------------------------------*/

	//同时计算原子的力、能量和维里
	void Non_Bond_14_LJ_CF_Force_With_Atom_Energy_And_Virial(const UINT_VECTOR_LJ_TYPE *uint_crd, const VECTOR scaler, const float *LJ_type_A, const float *LJ_type_B, VECTOR *frc, float *atom_energy, float *atom_virial);

	//获得能量
	float Get_14_LJ_Energy(const UINT_VECTOR_LJ_TYPE *uint_crd, const VECTOR scaler, const float *LJ_type_A, const float *LJ_type_B, int is_download = 1);
	float Get_14_CF_Energy(const UINT_VECTOR_LJ_TYPE *uint_crd, const VECTOR scaler, int is_download = 1);
	/*-----------------------------------------------------------------------------------------
	下面的函数是其他需求的排列组合，但是接口没有特地优化，如果自己需要，可能需要修改接口或写一个重载函数
	------------------------------------------------------------------------------------------*/

	void Non_Bond_14_LJ_Force(const UINT_VECTOR_LJ_TYPE *uint_crd, const VECTOR scaler, const float *LJ_type_A, const float *LJ_type_B, VECTOR *frc);
	void Non_Bond_14_LJ_CF_Force(const UINT_VECTOR_LJ_TYPE *uint_crd, const VECTOR scaler, const float *LJ_type_A, const float *LJ_type_B, VECTOR *frc);

	void Non_Bond_14_LJ_Energy(const UINT_VECTOR_LJ_TYPE *uint_crd, const VECTOR scaler, const float *LJ_type_A, const float *LJ_type_B);
	void Non_Bond_14_CF_Energy(const UINT_VECTOR_LJ_TYPE *uint_crd, const VECTOR scaler);
	void Non_Bond_14_LJ_CF_Energy(const UINT_VECTOR_LJ_TYPE *uint_crd, const VECTOR scaler, const float *LJ_type_A, const float *LJ_type_B);

	//能量列表导出并非给原子上，用于SITS等分能量方法
	void Non_Bond_14_LJ_Atom_Energy(const UINT_VECTOR_LJ_TYPE *uint_crd, const VECTOR scaler, const float *LJ_type_A, const float *LJ_type_B,float *atom_ene);
	void Non_Bond_14_CF_Atom_Energy(const UINT_VECTOR_LJ_TYPE *uint_crd, const VECTOR scaler,float *atom_ene);
	void Non_Bond_14_LJ_CF_Force_With_Atom_Energy(const UINT_VECTOR_LJ_TYPE *uint_crd, const VECTOR scaler, const float *LJ_type_A, const float *LJ_type_B, VECTOR *frc,float *atom_energy);

	
	
};

#endif //NB14_CUH(nb14.cuh)