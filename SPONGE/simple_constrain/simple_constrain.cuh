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
	float constrain_k;//���������˵�и�����ϵ�����̶������ǵ���ʱ���и�ϵ��k=m1*m2/(m1+m2)
};
struct SIMPLE_CONSTARIN
{
	char module_name[CHAR_LENGTH_MAX];
	int is_initialized = 0;
	int is_controller_printf_initialized = 0;
	int last_modify_date = 20210525;

	//Լ��������ʹ����ѭ���и��º��������ϸ��������ķ��������ǰ��pair����һ�£��������õ�����Լ�������ꡣ
	VECTOR *constrain_frc = NULL;
	//ÿ�Ե�ά��
	float *d_pair_virial = NULL;
	//��ά��
	float *d_virial = NULL;
	//����constrain���������еĲ���΢����ԭ��uint����
	UNSIGNED_INT_VECTOR *test_uint_crd = NULL;

	//��ѭ���и���ǰ��pair������Ϣ
	VECTOR *last_pair_dr = NULL;

	float dt_inverse;
	float half_exp_gamma_plus_half;//Ϊ0.5*(1.+exp_gamma)

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


	struct INFORMATION//����info�������Ϣ������Դ��md_info�������Initial_Simple_Constrainһ��ֱ�Ӵ��ν����õ��������������Դ
	{
		int atom_numbers = 0;
		float dt = 0.001f;
		VECTOR uint_dr_to_dr_cof;//��ϵ���ɽ��޷�����������֮���Ϊʵ������֮��
		VECTOR quarter_crd_to_uint_crd_cof;//��ϵ���ɽ�ʵ�������Ϊ��Ӧ��һ�볤�ȵ��޷�����������
		float volume; //���

		float exp_gamma = 1.0f;//��������ԡ��ʱ����Ҫ�޸����ֵ��Ϊ��ԡ��ʹ�õ�exp_gamma������֮�����1
		float step_length = 1.0f;//��������ʱѡȡ�Ĳ���������Ϊ1.���Ըպ��ϸ���������constrain
								//���������弰���ϵ����������Сһ�����ȶ�������֮��������Ҫ�������������
		int iteration_numbers = 25;//��������
		float constrain_mass = 3;//������С�ڸ�ֵ��ԭ�ӽ�������
	}info;

	//Ĭ�ϵ�Initial��Ҫ���������˳��
	//Add_HBond_To_Constrain_Pair
	//Add_HAngle_To_Constrain_Pair
	//Initial_Simple_Constrain

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
	//����exp_gammaΪ��֮��������ԡ��exp_gamma
	void Initial_Simple_Constrain
		(CONTROLLER *controller, const int atom_numbers, const float dt, const VECTOR box_length, const float exp_gamma, const int is_Minimization, float *atom_mass, int *system_freedom);
	
	//����ڴ�
	void Clear();
	//��¼����ǰ�ľ���
	void Remember_Last_Coordinates(VECTOR *crd);
	//����Լ������
	void Constrain
		(VECTOR *crd, VECTOR *vel, const float *mass_inverse, int need_pressure, float *d_pressure);
	//����仯ʱ�Ĳ�������
	void Update_Volume(VECTOR box_length);


};
#endif //SIMPLE_CONSTARIN_CUH(simple_constrain.cuh)
