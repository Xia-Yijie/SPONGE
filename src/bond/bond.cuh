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


#ifndef BOND_CUH
#define BOND_CUH
#include "../common.cuh"
#include "../control.cuh"

struct BOND
{
	char module_name[CHAR_LENGTH_MAX];
	int is_initialized = 0;
	int is_controller_printf_initialized = 0;
	int last_modify_date = 20210408;

	//E_bond=k*(|r_a-r_b|-r0)^2
	int bond_numbers = 0;
	int *h_atom_a = NULL;
	int *d_atom_a = NULL;
	int *h_atom_b = NULL;
	int *d_atom_b = NULL;
	float *d_k = NULL;
	float *h_k = NULL;
	float *d_r0 = NULL;
	float *h_r0 = NULL;

	float *h_bond_ene=NULL;
	float *d_bond_ene=NULL;
	float *d_sigma_of_bond_ene=NULL;
	float *h_sigma_of_bond_ene=NULL;

	//cuda计算分配相关参数
	int threads_per_block=128;

	//初始化模块
	void Initial(CONTROLLER *controller, char *module_name = NULL);

	//清空模块
	void Clear();
	//内存分配
	void Memory_Allocate();
	//从parm7文件中读取键信息 
	void Read_Information_From_AMBERFILE(const char *file_name, CONTROLLER controller);
	//拷贝cpu中的参数到gpu
	void Parameter_Host_To_Device();

	/*-----------------------------------------------------------------------------------------
	下面的函数是普通md的需求
	------------------------------------------------------------------------------------------*/

	//同时计算力，并将能量和维里加到每个原子头上
	void Bond_Force_With_Atom_Energy_And_Virial(const UNSIGNED_INT_VECTOR *unit_crd, const VECTOR scaler, VECTOR *frc, float *atom_energy, float *atom_virial);

	//获得能量
	float Get_Energy(const UNSIGNED_INT_VECTOR *unit_crd, const VECTOR scaler, int is_download = 1);

	/*-----------------------------------------------------------------------------------------
	下面的函数是其他需求的排列组合，但是接口没有特地优化，如果自己需要，可能需要修改接口或写一个重载函数
	------------------------------------------------------------------------------------------*/


	//将能量从gpu转移到cpu上
	void Energy_Device_To_Host();

	//计算bond force（整数坐标，缩放因子（将整数坐标变为实数坐标），原子受力）
	void Bond_Force(const UNSIGNED_INT_VECTOR *unit_crd, const VECTOR scaler, VECTOR *frc);

	//计算bond energy（整数坐标，缩放因子（将整数坐标变为实数坐标），原子受力）
	void Bond_Energy(const UNSIGNED_INT_VECTOR *unit_crd, const VECTOR scaler);

	//计算bond energy（整数坐标，缩放因子（将整数坐标变为实数坐标），原子能量）
	//将能量数组暴露出来，并将其加到每个原子的头上（用于SITS之类的分能过程）
	void Bond_Atom_Energy(const UNSIGNED_INT_VECTOR *unit_crd, const VECTOR scaler,float *atom_ene);
	
	//计算bond force并同时计算能量并加到原子能量列表上
	void Bond_Force_With_Atom_Energy(const UNSIGNED_INT_VECTOR *unit_crd, const VECTOR scaler, VECTOR *frc, float *atom_energy);
	
	//计算bond force并同时计算原子维里
	void Bond_Force_With_Atom_Virial(const UNSIGNED_INT_VECTOR *unit_crd, const VECTOR scaler, VECTOR *frc, float *atom_virial);


};

#endif //BOND (bond.cuh)
