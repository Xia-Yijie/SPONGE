#ifndef SOLVING_POISSON_EQUATION_CUH
#define SOLVING_POISSON_EQUATION_CUH
#include "../common.cuh"
struct FINITE_GRID_NEIGHBOR
{
	int x;//x增加方向
	int x_;//x减少方向
	int y;
	int y_;
	int z;
	int z_;
};

//两个列表交替SOR迭代，对于Nx Ny Nz全是偶数的情况（默认），该迭代得到的结果需要取平均才是最终解
__global__ void SOR_Iteration2(const int grid_numbers, const FINITE_GRID_NEIGHBOR *neighbor,
	const float *rho, const VECTOR scaler, const float gamma1, const float gamma2, float *phi, float *phi2);
//只对从start_grid_numbers到end_grid_numbers序号的格子进行差分迭代
__global__ void SOR_Iteration3(const int start_grid_numbers, const int end_grid_numbers, const FINITE_GRID_NEIGHBOR *neighbor,
	const float *rho, const VECTOR scaler, const float gamma1, const float gamma2, float *phi, float *phi2);

//求两个列表phi，phi2的平均并存入phi中
__global__ void SOR_Iteration2_Mean(const int grid_numbers, float *phi, float *phi2);

//计算前后两次phi值的距离，用于收敛判断
__global__ void Calculate_Phi_Convergence(const int element_numbers, const float *phi1, const float *phi2, float *sum);

#endif //SOLVING_POISSON_EQUATION_CUH