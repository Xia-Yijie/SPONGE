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


#ifndef MC_BAROSTAT_CUH
#define MC_BAROSTAT_CUH
#include "../common.cuh"
#include "../control.cuh"

//用于记录与计算MC控压相关的信息
struct MC_BAROSTAT_INFORMATION
{
	char module_name[CHAR_LENGTH_MAX];
	int is_initialized = 0;
	int is_controller_printf_initialized = 0;
	int last_modify_date = 20210827;

	int update_interval = 100;   //每多少步进行一次MC尝试
	int check_interval = 10;  //每多少次MC尝试以后进行一次DeltaV_max取值的检查（使得MC接受概率在40%~50%
	int reject = 1;  //接受与否
	VECTOR *frc_backup;     //备份力，以便还原
	VECTOR *crd_backup;     //备份坐标，以便还原
	int scale_coordinate_by_molecule = 1; //按照分子的质心进行规度坐标
	//成功率记录相关
	int total_count = 0;    //总共进行的MC尝试次数
	int accep_count = 0;    //接受的MC尝试次数
	float accept_rate = 0;  //接受率
	float accept_rate_low;  //接受率的低限
	float accept_rate_high; //接受率的高限

	
	//每次允许变化最大的体积
	float DeltaV_max = 200.0;

	//// 温度常数， 等于1 / k_B / T
	//float beta;

	//设定的压强
	float p0;

	//初始体积
	float V0;

	//每次尝试变化的体积
	double DeltaV;
	double newV;
	double VDevided;

	//坐标系数因子
	double crd_scale_factor;

	//概率项
	float energy_old;
	float energy_new;
	float extra_term;
	float accept_possibility;

	//初始化
	void Initial(CONTROLLER *controller, int atom_numbers,
		float target_pressure, VECTOR boxlength, int res_is_initialized, char *module_name = NULL);
	
	//输入盒子信息，得到坐标变化因子
	void Volume_Change_Attempt(VECTOR boxlength);
	
	//判断是否接受
	int Check_MC_Barostat_Accept();

	//对猜测值根据接受率更新
    void Delta_V_Max_Update();

	//判断步数条件，决定是否让计算力的时候也计算能量
	void Ask_For_Calculate_Potential(int steps, int *need_potential);
	
};

#endif //MC_BAROSTAT_CUH(MC_barostat.cuh)
