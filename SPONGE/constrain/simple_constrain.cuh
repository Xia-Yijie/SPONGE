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
#include "constrain.cuh"

struct SIMPLE_CONSTRAIN
{
	char module_name[CHAR_LENGTH_MAX];
	int is_initialized = 0;
	int is_controller_printf_initialized = 0;
	int last_modify_date = 20211222;

	CONSTRAIN *constrain;

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

	float step_length = 1.0f;//��������ʱѡȡ�Ĳ���������Ϊ1.���Ըպ��ϸ���������constrain
							//���������弰���ϵ����������Сһ�����ȶ�������֮��������Ҫ�������������
	int iteration_numbers = 25;//��������
	
	//�ڼ������constrain_pair���ʼ��
	//����exp_gammaΪ��֮��������ԡ��exp_gamma
	void Initial_Simple_Constrain
		(CONTROLLER *controller, CONSTRAIN *constrain, const char *module_name = NULL);
	
	//����ڴ�
	void Clear();
	//��¼����ǰ�ľ���
	void Remember_Last_Coordinates(VECTOR *crd, UNSIGNED_INT_VECTOR* uint_crd, VECTOR scaler);
	//����Լ������
	void Constrain
		(VECTOR *crd, VECTOR *vel, const float *mass_inverse, const float *d_mass, VECTOR box_length, int need_pressure, float *d_pressure);
	//����仯ʱ�Ĳ�������
	void Update_Volume(VECTOR box_length);

};



#endif //SIMPLE_CONSTARIN_CUH(simple_constrain.cuh)
