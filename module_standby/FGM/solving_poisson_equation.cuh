#ifndef SOLVING_POISSON_EQUATION_CUH
#define SOLVING_POISSON_EQUATION_CUH
#include "../common.cuh"
struct FINITE_GRID_NEIGHBOR
{
	int x;//x���ӷ���
	int x_;//x���ٷ���
	int y;
	int y_;
	int z;
	int z_;
};

//�����б���SOR����������Nx Ny Nzȫ��ż���������Ĭ�ϣ����õ����õ��Ľ����Ҫȡƽ���������ս�
__global__ void SOR_Iteration2(const int grid_numbers, const FINITE_GRID_NEIGHBOR *neighbor,
	const float *rho, const VECTOR scaler, const float gamma1, const float gamma2, float *phi, float *phi2);
//ֻ�Դ�start_grid_numbers��end_grid_numbers��ŵĸ��ӽ��в�ֵ���
__global__ void SOR_Iteration3(const int start_grid_numbers, const int end_grid_numbers, const FINITE_GRID_NEIGHBOR *neighbor,
	const float *rho, const VECTOR scaler, const float gamma1, const float gamma2, float *phi, float *phi2);

//�������б�phi��phi2��ƽ��������phi��
__global__ void SOR_Iteration2_Mean(const int grid_numbers, float *phi, float *phi2);

//����ǰ������phiֵ�ľ��룬���������ж�
__global__ void Calculate_Phi_Convergence(const int element_numbers, const float *phi1, const float *phi2, float *sum);

#endif //SOLVING_POISSON_EQUATION_CUH