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


#ifndef COORDINATE_MOLECULAR_MAP
#define COORDINATE_MOLECULAR_MAP
#include "../common.cuh"
#include "../control.cuh"
#include <deque>

//20210420 ����ӳ��Ŀǰ�����ų������ġ����������ų����ԭ�ӻᱻ��Ϊ����ͬһ������
struct CoordinateMolecularMap
{
	char module_name[CHAR_LENGTH_MAX];
	int is_initialized = 0;
	int last_modify_date = 20210830;

	//��ϵ������Ϣ
	int atom_numbers=0;
	VECTOR box_length;

	//������Ҫ����������Ľ�wrap�������¼
	VECTOR *nowrap_crd = NULL;
	VECTOR *old_crd = NULL;
	INT_VECTOR *box_map_times = NULL;
	VECTOR *h_nowrap_crd = NULL;
	VECTOR *h_old_crd = NULL;
	INT_VECTOR *h_box_map_times = NULL;

	int threads_per_block = 256;
	int blocks_per_grid = 20;
	//ע�⴫���crd��device�ϵ�ַ��һ���ʼ����ʱ�����������������
	void Initial(int atom_numbers, VECTOR box_length, VECTOR *crd, 
		const int exclude_numbers, const int *exclude_length, const int *exclude_start, const int *exclude_list, char *module_name = NULL);
	//����ڴ�
	void Clear();

	//����ʵʱģ���е�crd���꣨GPU�ϣ����������nowarp_crd�����ڷ������������
	void Calculate_No_Wrap_Crd(const VECTOR *crd);
	//��ÿ���ж�ԭ�ӽ���������ӳ��Ĳ��������������������Ը��¿���Ӵ�����¼��
	void Refresh_BoxMapTimes(const VECTOR *crd);

	//CPU�Ϻ������жϽ��ڱ�ģ���н���box_map����ǰ���old_crd,crd�Ĵ�Խ������Ϣbox_map_times
	void Record_Box_Map_Times_Host(int atom_numbers, VECTOR *crd, VECTOR *old_crd, INT_VECTOR *box_map_times, VECTOR box);

	void Update_Volume(VECTOR box_length);

};
#endif //COORDINATE_MOLECULAR_MAP(crd_molecular_map.cuh)