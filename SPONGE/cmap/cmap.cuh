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


#ifndef CMAP_CUH
#define CMAP_CUH
#include "../common.cuh"
#include "../control.cuh"


struct CMAP
{
	char module_name[CHAR_LENGTH_MAX];
	int is_initialized = 0;
	int is_controller_printf_initialized = 0;
	int last_modify_date = 20211129;

	//�����ļ���Ϣ����
	int tot_cmap_num = 0;
	int uniq_cmap_num = 0;
	int tot_gridpoint_num = 0;
	int uniq_gridpoint_num = 0;
	int* cmap_type = NULL;
	int* cmap_resolution = NULL;
	float* grid_value = NULL;

	//cuda���������ز���
	int threads_per_block = 128;

	//��ֵϵ�������������൱�ڽ����Է�����õ���ֵ����ʽϵ��
	/*
	ÿ16����������һ�β�ֵ����Ԫ����ʽ��ʽΪ��F(x,y) = \sum_{i.j=0,1,2,3}(c_{ij}x^iy^j),����ȡ��Ϊ4*4,����PBC��
	�Գ�ʼ���ݸ�����߳��ȹ�һ���󣬿��Եõ�ϵ������(16*1):

	                     c = A^{-1}p
	
	����p�����ǰ������ֵ(4)��һ�ײ��(8)�����ײ��(4)�ģ�16*1������
	*/

	//���ӵ������˫�������������ļ�

	int *h_atom_a = NULL;
	int *d_atom_a = NULL;
	int *h_atom_b = NULL;
	int *d_atom_b = NULL;
	int *h_atom_c = NULL;
	int *d_atom_c = NULL;
	int *h_atom_d = NULL;
	int *d_atom_d = NULL;
	int *h_atom_e = NULL;
	int *d_atom_e = NULL;
	
	int* d_cmap_resolution = NULL;
	int* d_cmap_type = NULL;
	float* inter_coeff = NULL;
	float *d_inter_coeff = NULL;
	float *d_cmap_ene = NULL;
	float *h_sigma_of_cmap_ene = NULL;
	float *d_sigma_of_cmap_ene = NULL;

	float *h_cmap_force = NULL;
	float *d_cmap_force = NULL;//cmap_test_temp



	//��ʼ��ģ��
	void Initial(CONTROLLER *controller, const char *module_name = NULL);
	//���ģ��
	void Clear();
	//�ڴ����
	void Memory_Allocate();
	//��parm7�ļ��ж�ȡ����Ϣ 
	void Read_Information_From_AMBERFILE(const char *file_name, CONTROLLER controller);
	//�����ֵϵ��
	void Interpolation(int* resolution, float *grid_value,CONTROLLER controller);
	//CUDA����
	void Parameter_Host_to_Device();

	

	//������������
	void CMAP_Force_with_Atom_Energy(const UNSIGNED_INT_VECTOR *uint_crd, const VECTOR scaler, VECTOR *frc, float *atom_energy);
	float Get_Energy(const UNSIGNED_INT_VECTOR *uint_crd, const VECTOR scaler, int is_download = 1);



};
#endif
























