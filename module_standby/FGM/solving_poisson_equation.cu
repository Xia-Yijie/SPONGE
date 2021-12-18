#include "solving_poisson_equation.cuh"
__global__ void SOR_Iteration2(const int grid_numbers, const FINITE_GRID_NEIGHBOR *neighbor,
	const float *rho, const VECTOR scaler, const float gamma1, const float gamma2, float *phi, float *phi2)
{
	int grid_i = blockDim.x*blockIdx.x + threadIdx.x;
	if (grid_i < grid_numbers)
	{

		FINITE_GRID_NEIGHBOR lin = neighbor[grid_i];
		float phi_lin = scaler.x*(phi[lin.x] + phi[lin.x_]);
		phi_lin = phi_lin + scaler.y*(phi[lin.y] + phi[lin.y_]);
		phi_lin = phi_lin + scaler.z*(phi[lin.z] + phi[lin.z_]);

		phi_lin = phi_lin + rho[grid_i];
		phi2[grid_i] = gamma1*phi2[grid_i] + gamma2*phi_lin;

	}
}
__global__ void SOR_Iteration3(const int start_grid_numbers, const int end_grid_numbers, const FINITE_GRID_NEIGHBOR *neighbor,
	const float *rho, const VECTOR scaler, const float gamma1, const float gamma2, float *phi, float *phi2)
{
	int grid_i = blockDim.x*blockIdx.x + threadIdx.x + start_grid_numbers;
	if (grid_i < end_grid_numbers)
	{

		FINITE_GRID_NEIGHBOR lin = neighbor[grid_i];
		float phi_lin = scaler.x*(phi[lin.x] + phi[lin.x_]);
		phi_lin = phi_lin + scaler.y*(phi[lin.y] + phi[lin.y_]);
		phi_lin = phi_lin + scaler.z*(phi[lin.z] + phi[lin.z_]);

		phi_lin = phi_lin + rho[grid_i];
		phi2[grid_i] = gamma1*phi2[grid_i] + gamma2*phi_lin;

	}
}

__global__ void SOR_Iteration2_Mean(const int grid_numbers, float *phi, float *phi2)
{
	int grid_i = (blockDim.x*blockIdx.x + threadIdx.x);
	if (grid_i < grid_numbers)
	{
		phi[grid_i] = 0.5*(phi[grid_i] + phi2[grid_i]);
		phi2[grid_i] = phi[grid_i];
	}
}

__global__ void Calculate_Phi_Convergence(const int element_numbers, const float *phi1, const float *phi2, float *sum)
{
	if (threadIdx.x == 0)
	{
		sum[0] = 0;
	}
	__syncthreads();
	float lin = 0;
	for (int i = threadIdx.x; i < element_numbers; i = i + blockDim.x)
	{
		lin = lin + (phi1[i] - phi2[i])*(phi1[i] - phi2[i]);
	}
	atomicAdd(sum, lin);
}