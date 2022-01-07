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


#ifndef CONSTARIN_CUH
#define CONSTARIN_CUH
#include "../common.cuh"
#include "../control.cuh"

struct CONSTRAIN_PAIR
{
	int atom_i_serial;
	int atom_j_serial;
	float constant_r;
	float constrain_k;//���������˵�и�����ϵ�����̶������ǵ���ʱ���и�ϵ��k=m1*m2/(m1+m2)
};

struct CONSTRAIN
{
	char module_name[CHAR_LENGTH_MAX];
	int is_initialized = 0;
	int is_controller_printf_initialized = 0;
	int last_modify_date = 20211222;

	int atom_numbers = 0;
	float dt = 0.001f;
	float dt_inverse;
	VECTOR uint_dr_to_dr_cof;//��ϵ���ɽ��޷�����������֮���Ϊʵ������֮��
	VECTOR quarter_crd_to_uint_crd_cof;//��ϵ���ɽ�ʵ�������Ϊ��Ӧ��һ�볤�ȵ��޷�����������
	float volume; //���

	float v_factor = 1.0f;  //һ�����ֲ���,һ��΢С����F���ٶȵ�Ӱ�죬��dv = v_factor * F * dt/m
	float x_factor = 1.0f;  //һ�����ֲ���,һ��΢С����F��λ�Ƶ�Ӱ�죬��dx = x_factor * F * dt * dt/m 
	float constrain_mass = 3.3;//������С�ڸ�ֵ��ԭ�ӽ�������

	//�ڳ�ʼ����ʱ���õ�����ʵ�ʼ����в���ʹ��,�ڳ�ʼ��ʱ�Ѿ����ͷ�
	int bond_constrain_pair_numbers = 0;
	int angle_constrain_pair_numbers = 0;
	CONSTRAIN_PAIR *h_bond_pair = NULL;
	CONSTRAIN_PAIR *h_angle_pair = NULL;

	//��ʵ�ʼ�����ʹ�ã���ϵ�ܵ�constrain pair
	int constrain_pair_numbers = 0;
	CONSTRAIN_PAIR *constrain_pair = NULL;
	CONSTRAIN_PAIR *h_constrain_pair = NULL;

	//������ʱ��¼bond����Ϣ������angle������bond����
	//��Щָ��ָ��Ŀռ䲢���ɱ�ģ�������Ҳ��ɱ�ģ���ͷ�
	struct BOND_INFORMATION
	{
		int bond_numbers;
		const int *atom_a = NULL;
		const int *atom_b = NULL;
		const float *bond_r = NULL;
	}bond_info;


	//Ĭ�ϵ�Initial��Ҫ���������˳��
	//Add_HBond_To_Constrain_Pair
	//Add_HAngle_To_Constrain_Pair
	//Initial_Constrain

	//20201125 ����MD_INFORMATION������ʱû�м���ԭ�����������Բ���ԭ�������ж�H������ֱ�ӱȽ�����
	//��������Сʱ������Ϊ��H
	//�����ָ��ָ��HOST�ڴ�
	void Add_HBond_To_Constrain_Pair
		(CONTROLLER *controller, const int bond_numbers, const int *atom_a, const int *atom_b, const float *bond_r,
		const float *atom_mass, char *module_name = NULL);//Ҫ�����ָ��host���ڴ��ָ��

	//��Ҫ��������Add_HBond_To_Constrain_Pair֮��������
	//�����ָ��ָ��HOST�ڴ�
	void Add_HAngle_To_Constrain_Pair
		(CONTROLLER *controller, const int angle_numbers, const int *atom_a, const int *atom_b, const int *atom_c,
		const float *angle_theta, const float *atom_mass);//Ҫ�����ָ��host���ڴ��ָ��
	
	//�ڼ������constrain_pair���ʼ��
	//�м��exp_gammaΪ��֮��������ԡ��exp_gamma
	void Initial_Constrain
		(CONTROLLER *controller, const int atom_numbers, const float dt, const VECTOR box_length, const float exp_gamma, const int is_Minimization, float *atom_mass, int *system_freedom);
	
	//����ڴ�
	void Clear();
	//��¼����ǰ�ľ���

	void Update_Volume(VECTOR box_length);

};



#endif //SIMPLE_CONSTARIN_CUH(simple_constrain.cuh)
