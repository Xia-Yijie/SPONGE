#ifndef FGM_COMMON_CUH
#define FGM_COMMON_CUH
#include "solving_poisson_equation.cuh"
#include "cufft.h"

enum FGM_MODE
{
	BULK_MODE = 0x00000001,
	INFINITY_Z_MODE = 0x00000002,

	READ_POINT_PHI_DAT = 0x00000010,
	CALCULATE_POINT_PHI_DAT = 0x00000020,

	SAME_FINITE_DIMENSION = 0x00000100,
	DIFFERENT_FINITE_DIMENSION = 0x00000200,

	GREEN_SAMPLING_TEXTURE_VERSION=0X00001000
};
struct FINITE_GREEN_METHOD
{
private:
	CONSTANT const_parameter;

	int grid_numbers=0;
	float grid_numbers_inverse;
	int layer_numbers = 0;
	INT_VECTOR grid_dimension ;
	VECTOR grid_length;
	VECTOR uint_crd_to_grid_serial;
	FINITE_GRID_NEIGHBOR *h_neighbor = NULL;
	FINITE_GRID_NEIGHBOR *neighbor = NULL;
	int Initial_Grid_Neighbor();

	float *point_phi = NULL;//点电荷生成的电势分布
	float *h_point_phi = NULL;//点电荷生成的电势分布
	cufftComplex *fft_of_point_phi = NULL;
	float *phi = NULL;//体系的电势分布
	float *h_phi = NULL;//体系的电势分布
	cufftComplex *fft_of_phi = NULL;
	float *rho = NULL;//体系的电荷分布
	cufftComplex *fft_of_rho = NULL;

	cufftHandle plan_r2c;
	cufftHandle plan_c2r;

	//texture 内存加速的格林积分点采样相关
	VECTOR uint_crd_scale_to_1;
	VECTOR half_grid_lenth_scale_to_1;
	cudaArray* cuArray_phi=NULL;
	cudaMemcpy3DParms copyParams_phi;
	cudaTextureObject_t texObj_phi;
	int Initial_TextureObj_Phi();

	struct GREEN_SPHERE_INFORMATION
	{
		int point_numbers=0.;//格林球面离散点数目
		float sphere_radius=10.;//格林球面半径，一般与分子模拟的cutoff相同
		float ds=0.;//格林球面积分前面的系数，也可认为是离散的面元大小，大小为4pi/3/sphere_radius^2/point_numbers
		VECTOR *h_crd=NULL;//离散点的真实三维坐标，球面半径为1.
		VECTOR *crd=NULL;//离散点的真实三维坐标,球面半径为sphere_radius
		VECTOR *crd_with_ds=NULL;//已经乘面元大小的球面离散点坐标
		INT_VECTOR *int_crd=NULL;//离散点的整数坐标（存在负数，由周期边界大小决定）
	}gs_info;
	struct FINITE_DIFFERENCE_SOLVING_POISSON_EQUATION_INFORMATION
	{
		int grid_level = 1;
		int fd_grid_numbers = 0;//用有限差分法求解点电势分布的时候使用的格点数目（只允许为grid_numbersgrid_level3次方倍）
		float first_gamma = 1.;//SOR超松弛迭代的gamma因子，用于全局近似收敛
		int first_iteration_steps = 0;//全局近似迭代次数
		float second_gamma = 1.;//SOR超松弛迭代的gamma因子，用于全局近似收敛后的精确收敛
		int second_iteration_steps = 0;//精确收敛迭代次数

		VECTOR fd_grid_length;
		VECTOR fd_grid_length_inverse_square;
		FINITE_GRID_NEIGHBOR *fd_neighbor = NULL;
		FINITE_GRID_NEIGHBOR *h_fd_neighbor = NULL;
		INT_VECTOR fd_grid_dimension;
		int fd_layer_numbers;
		float SOR_gamma = 1.;
		float SOR_gamma_1 = 0.;
		float SOR_gamma_2 = 1.;
		int Refresh_Gamma(float target_gamma);
		int Solving_Poisson_Equation(float *h_point_phi);
		int Initial_Grid_Neighbor();
	}fd_info;

public:
	struct FINITE_GREEN_METHOD_CONTROL_INFORMATION
	{
		int fgm_mode = 0;
		VECTOR box_length;
		VECTOR target_grid_length;//预期FGMD求解泊松方程所划分格子的尺寸
		char sphere_pos_file_name[128];//记录采样点坐标的文件
		char point_phi_file_name[128];//记录点电荷电势分布的文件

		int sphere_point_numbers=0;
	}fgm_info;
	int Set_Finite_Difference_Parameter(int grid_level,float first_gamma,int first_iteration_steps,float second_gamma,int second_iteration_steps);
	int Initial_FGM();
	int Calculate_FGM_Force
		(const int charge_numbers,
		const float *charge_density_list,
		const UNSIGNED_INT_VECTOR *uint_crd,
		VECTOR *frc);
	int Export_Phi(char *file_name);
};
#endif //FGM_COMMON_CUH