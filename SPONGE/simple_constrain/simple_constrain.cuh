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


#ifndef SIMPLE_CONSTARIN_CUH
#define SIMPLE_CONSTARIN_CUH
#include "../common.cuh"
#include "../control.cuh"

struct CONSTRAIN_PAIR
{
	int atom_i_serial;
	int atom_j_serial;
	float constant_r;
	float constrain_k;//这个并不是说有个弹性系数来固定，而是迭代时，有个系数k=m1*m2/(m1+m2)
};
struct SIMPLE_CONSTARIN
{
	char module_name[CHAR_LENGTH_MAX];
	int is_initialized = 0;
	int is_controller_printf_initialized = 0;
	int last_modify_date = 20210525;

	//约束内力，使得主循环中更新后的坐标加上该力（力的方向与更新前的pair方向一致）修正，得到满足约束的坐标。
	VECTOR *constrain_frc = NULL;
	//每对的维里
	float *d_pair_virial = NULL;
	//总维里
	float *d_virial = NULL;
	//进行constrain迭代过程中的不断微调的原子uint坐标
	UNSIGNED_INT_VECTOR *test_uint_crd = NULL;

	//主循环中更新前的pair向量信息
	VECTOR *last_pair_dr = NULL;

	float dt_inverse;
	float half_exp_gamma_plus_half;//为0.5*(1.+exp_gamma)

	//在初始化的时候用到，在实际计算中不会使用,在初始化时已经被释放
	int bond_constrain_pair_numbers = 0;
	int angle_constrain_pair_numbers = 0;
	CONSTRAIN_PAIR *h_bond_pair = NULL;
	CONSTRAIN_PAIR *h_angle_pair = NULL;

	//在实际计算中使用，体系总的constrain pair
	int constrain_pair_numbers = 0;
	CONSTRAIN_PAIR *constrain_pair = NULL;
	CONSTRAIN_PAIR *h_constrain_pair = NULL;

	//用于暂时记录bond的信息，便于angle中搜索bond长度
	//这些指针指向的空间并不由本模块申请且不由本模块释放
	struct BOND_INFORMATION
	{
		int bond_numbers;
		const int *atom_a = NULL;
		const int *atom_b = NULL;
		const float *bond_r = NULL;
	}bond_info;


	struct INFORMATION//由于info里面的信息几乎来源于md_info，因此在Initial_Simple_Constrain一步直接传参进来得到，而无需额外来源
	{
		int atom_numbers = 0;
		float dt = 0.001f;
		VECTOR uint_dr_to_dr_cof;//该系数可将无符号整型坐标之差变为实际坐标之差
		VECTOR quarter_crd_to_uint_crd_cof;//该系数可将实际坐标变为对应的一半长度的无符号整型坐标
		float volume; //体积

		float exp_gamma = 1.0f;//用刘剑热浴的时候，需要修改这个值变为热浴中使用的exp_gamma，除此之外就是1
		float step_length = 1.0f;//迭代求力时选取的步长，步长为1.可以刚好严格求得两体的constrain
								//但对于三体及以上的情况，步长小一点会更稳定，但随之而来可能要求迭代次数增加
		int iteration_numbers = 25;//迭代步数
		float constrain_mass = 3;//对质量小于该值的原子进行限制
	}info;

	//默认的Initial需要按照下面的顺序：
	//Add_HBond_To_Constrain_Pair
	//Add_HAngle_To_Constrain_Pair
	//Initial_Simple_Constrain

	//20201125 由于MD_INFORMATION里面暂时没有加入原子序数，所以不用原子序数判断H，而是直接比较质量
	//当质量较小时，就认为是H
	//传入的指针指向HOST内存
	void Add_HBond_To_Constrain_Pair
		(CONTROLLER *controller, const int bond_numbers, const int *atom_a, const int *atom_b, const float *bond_r,
		const float *atom_mass, char *module_name = NULL);//要求均是指向host上内存的指针

	//需要在先运行Add_HBond_To_Constrain_Pair之后再运行
	//传入的指针指向HOST内存
	void Add_HAngle_To_Constrain_Pair
		(CONTROLLER *controller, const int angle_numbers, const int *atom_a, const int *atom_b, const int *atom_c,
		const float *angle_theta, const float *atom_mass);//要求均是指向host上内存的指针
	
	//在加入各种constrain_pair后初始化
	//最后的exp_gamma为朗之万刘剑热浴的exp_gamma
	void Initial_Simple_Constrain
		(CONTROLLER *controller, const int atom_numbers, const float dt, const VECTOR box_length, const float exp_gamma, const int is_Minimization, float *atom_mass, int *system_freedom);
	
	//清除内存
	void Clear();
	//记录更新前的距离
	void Remember_Last_Coordinates(VECTOR *crd);
	//进行约束迭代
	void Constrain
		(VECTOR *crd, VECTOR *vel, const float *mass_inverse, int need_pressure, float *d_pressure);
	//体积变化时的参数更新
	void Update_Volume(VECTOR box_length);


};
#endif //SIMPLE_CONSTARIN_CUH(simple_constrain.cuh)
