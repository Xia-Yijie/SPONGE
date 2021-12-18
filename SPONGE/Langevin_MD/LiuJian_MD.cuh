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


#ifndef LIUJIAN_MD_CUH
#define LIUJIAN_MD_CUH
#include "../common.cuh"
#include "../control.cuh"

struct LIUJIAN_MD_INFORMATION
{
	char module_name[CHAR_LENGTH_MAX];
	int is_initialized = 0;
	int is_controller_printf_initialized = 0;
	int last_modify_date = 20210525;

	int threads_per_block = 32;

	int atom_numbers;
	float dt;
	float half_dt;
	float gamma_ln;//碰撞频率
	float target_temperature;//热浴温度

	float exp_gamma;//刘剑动力学中参数( = expf(-gamma_ln*dt));
	int float4_numbers;//存储随机数的长度
	curandStatePhilox4_32_10_t *rand_state = NULL;//用于记录随机数发生器状态
	VECTOR *random_force = NULL;//存储随机力矢量，要求该数组的长度要能整除4且大于等于atom_numbers
	float *d_sqrt_mass = NULL;//用于刘剑热浴过程中随机力的原子等效质量
	float *h_sqrt_mass = NULL;//用于刘剑热浴过程中随机力的原子等效质量
	float *d_mass_inverse = NULL;//质量的倒数

	//使用速度上限的迭代方法而非盲目加大摩擦、降低温度、减少步长
	float max_velocity;

	//初始化（质量直接读取文件）
	void Initial(CONTROLLER *controller, float target_temperature, const char *module_name = NULL);
	//从AMBER中读取质量
	void Initial_From_AMBER(float **h_mass, const char *file, CONTROLLER controller);
	//清除内存
	void Clear();
	//迭代算法
	void MD_Iteration_Leap_Frog(VECTOR *frc, VECTOR *vel, VECTOR *acc, VECTOR *crd);

};


#endif //LIUJIAN_MD_CUH(LiuJian_MD.cuh)
