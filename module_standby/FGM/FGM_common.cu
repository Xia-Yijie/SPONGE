#include "FGM_common.cuh"
//一些辅助函数
int Check_2357_Factor(int number);
int Get_FFT_Patameter(float length, float factor);
__global__ void cufft_times(int element_numbers, cufftComplex *a, cufftComplex *b, cufftComplex *c);
__global__ void Distribute_Charge_To_Vertex
(const int atom_numbers, const float *charge_density, const UNSIGNED_INT_VECTOR *uint_crd,
const VECTOR scaler, const INT_VECTOR Nmax, const int Nxy,
const FINITE_GRID_NEIGHBOR *neighbor, float *rho);
__global__ void Green_Force_Texture_Ver
(const int atom_numbers, const UNSIGNED_INT_VECTOR *uint_crd, const float *atom_charge,
const int point_numbers, const INT_VECTOR *sphere_point, const VECTOR *sphere_direct, const VECTOR scaler,
const VECTOR half_grid_lenth_devide_box_length,
const cudaTextureObject_t phi, VECTOR *frc)
{
	int atom_i = blockDim.x*blockIdx.x + threadIdx.x;
	if (atom_i < atom_numbers)
	{
		float charge = atom_charge[atom_i];
		UNSIGNED_INT_VECTOR crd_i = uint_crd[atom_i];
		UNSIGNED_INT_VECTOR point_j;
		VECTOR r;
		VECTOR frc_lin = { 0., 0., 0. };
		float phi_lin;

		for (int j = 0; j < point_numbers; j = j + 1)
		{
			point_j.uint_x = crd_i.uint_x + sphere_point[j].int_x;
			point_j.uint_y = crd_i.uint_y + sphere_point[j].int_y;
			point_j.uint_z = crd_i.uint_z + sphere_point[j].int_z;
			r.x = (float)scaler.x*point_j.uint_x;
			r.y = (float)scaler.y*point_j.uint_y;
			r.z = (float)scaler.z*point_j.uint_z;//0~1 坐标

			r.x = r.x + half_grid_lenth_devide_box_length.x;
			r.y = r.y + half_grid_lenth_devide_box_length.y;
			r.z = r.z + half_grid_lenth_devide_box_length.z;

			phi_lin = tex3D<float>(phi, r.x, r.y, r.z);

			frc_lin.x = frc_lin.x - phi_lin*sphere_direct[j].x;
			frc_lin.y = frc_lin.y - phi_lin*sphere_direct[j].y;
			frc_lin.z = frc_lin.z - phi_lin*sphere_direct[j].z;

		}//green sphere j cycle
		frc_lin.x = frc_lin.x*charge;
		frc_lin.y = frc_lin.y*charge;
		frc_lin.z = frc_lin.z*charge;
		atomicAdd(&frc[atom_i].x, frc_lin.x);
		atomicAdd(&frc[atom_i].y, frc_lin.y);
		atomicAdd(&frc[atom_i].z, frc_lin.z);
	}
}
int FINITE_GREEN_METHOD::Initial_TextureObj_Phi()
{
	uint_crd_scale_to_1.x = (float)uint_crd_to_grid_serial.x / grid_dimension.int_x;
	uint_crd_scale_to_1.y = (float)uint_crd_to_grid_serial.y / grid_dimension.int_y;
	uint_crd_scale_to_1.z = (float)uint_crd_to_grid_serial.z / grid_dimension.int_z;

	half_grid_lenth_scale_to_1.x = 0.5*grid_length.x / fgm_info.box_length.x;
	half_grid_lenth_scale_to_1.y = 0.5*grid_length.y / fgm_info.box_length.y;
	half_grid_lenth_scale_to_1.z = 0.5*grid_length.z / fgm_info.box_length.z;

	// Allocate CUDA array in device memory
	cudaChannelFormatDesc channelDesc =
		cudaCreateChannelDesc(32, 0, 0, 0,
		cudaChannelFormatKindFloat);
	cudaExtent cuEx;
	cuEx.depth = grid_dimension.int_z;
	cuEx.height = grid_dimension.int_y;
	cuEx.width = grid_dimension.int_x;
	cudaMalloc3DArray(&cuArray_phi, &channelDesc, cuEx);
	// copy data to 3D array
	copyParams_phi = { 0 };
	copyParams_phi.srcPtr = make_cudaPitchedPtr((void *)phi, cuEx.width*sizeof(float), cuEx.width, cuEx.height);
	copyParams_phi.dstArray = cuArray_phi;
	copyParams_phi.extent = cuEx;
	copyParams_phi.kind = cudaMemcpyDeviceToDevice;
	cudaMemcpy3D(&copyParams_phi);

	// Specify texture
	struct cudaResourceDesc resDesc;
	memset(&resDesc, 0, sizeof(resDesc));
	resDesc.resType = cudaResourceTypeArray;
	resDesc.res.array.array = cuArray_phi;

	// Specify texture object parameters
	struct cudaTextureDesc texDesc;
	memset(&texDesc, 0, sizeof(texDesc));
	texDesc.addressMode[0] = cudaAddressModeWrap;
	texDesc.addressMode[1] = cudaAddressModeWrap;
	texDesc.addressMode[2] = cudaAddressModeWrap;
	texDesc.filterMode = cudaFilterModeLinear;
	texDesc.readMode = cudaReadModeElementType;
	texDesc.normalizedCoords = 1;

	// Create texture object
	cudaCreateTextureObject(&texObj_phi, &resDesc, &texDesc, NULL);
	return 1;
}
int FINITE_GREEN_METHOD::Initial_FGM()
{
	//体相模式
	if ((fgm_info.fgm_mode&0x0000000F)==BULK_MODE)
	{
		printf("Initial Finite Green Method:\n");
		//用给定的目标格点间距生成离散格点
		grid_dimension.int_x = (int)Get_FFT_Patameter(fgm_info.box_length.x, fgm_info.target_grid_length.x);
		grid_dimension.int_y = (int)Get_FFT_Patameter(fgm_info.box_length.y, fgm_info.target_grid_length.y);
		grid_dimension.int_z = (int)Get_FFT_Patameter(fgm_info.box_length.z, fgm_info.target_grid_length.z);
		layer_numbers = grid_dimension.int_x*grid_dimension.int_y;
		grid_numbers = grid_dimension.int_z*layer_numbers;
		grid_numbers_inverse =(float) 1. / grid_numbers;
		grid_length.x = (float)fgm_info.box_length.x / grid_dimension.int_x;
		grid_length.y = (float)fgm_info.box_length.y / grid_dimension.int_y;
		grid_length.z = (float)fgm_info.box_length.z / grid_dimension.int_z;
		printf("	GridDimension:%8d,%8d,%8d\n", grid_dimension.int_x, grid_dimension.int_y, grid_dimension.int_z);
		printf("	GridLength:%8f,%8f,%8f\n", grid_length.x, grid_length.y, grid_length.z);
		Cuda_Malloc_Safely((void**)&fft_of_point_phi, sizeof(cufftComplex)* grid_numbers);
		Cuda_Malloc_Safely((void**)&fft_of_phi, sizeof(cufftComplex)* grid_numbers);
		Cuda_Malloc_Safely((void**)&fft_of_rho, sizeof(cufftComplex)* grid_numbers);
		Cuda_Malloc_Safely((void**)&point_phi, sizeof(float)* grid_numbers);
		Cuda_Malloc_Safely((void**)&phi, sizeof(float)* grid_numbers);
		Malloc_Safely((void**)&h_phi, sizeof(float)* grid_numbers);
		Cuda_Malloc_Safely((void**)&rho, sizeof(float)* grid_numbers);

		uint_crd_to_grid_serial.x = fgm_info.box_length.x / const_parameter.uint_max_float / grid_length.x;
		uint_crd_to_grid_serial.y = fgm_info.box_length.y / const_parameter.uint_max_float / grid_length.y;
		uint_crd_to_grid_serial.z = fgm_info.box_length.z / const_parameter.uint_max_float / grid_length.z;

		Malloc_Safely((void**)&h_point_phi, sizeof(float)* grid_numbers);

		cufftPlan3d(&plan_r2c, grid_dimension.int_x, grid_dimension.int_y, grid_dimension.int_z, CUFFT_R2C);
		cufftPlan3d(&plan_c2r, grid_dimension.int_x, grid_dimension.int_y, grid_dimension.int_z, CUFFT_C2R);
		
		Initial_Grid_Neighbor();

		//从已经用有限元计算好的点电荷电势分布文件中获得点电荷电势分布
		if ((fgm_info.fgm_mode & 0x000000F0) == READ_POINT_PHI_DAT)
		{
			//有限元计算点电荷分布的时候采用与分子模拟做FFT时同样的格点数目
			if ((fgm_info.fgm_mode & 0x00000F00) == SAME_FINITE_DIMENSION)
			{
				FILE *point_phi_in = NULL;
				Open_File_Safely(&point_phi_in, fgm_info.point_phi_file_name, "rb");
				fread(h_point_phi, sizeof(float), grid_numbers, point_phi_in);
				cudaMemcpy(point_phi, h_point_phi, sizeof(float)*grid_numbers, cudaMemcpyHostToDevice);
				fclose(point_phi_in);
				printf("	Read Point_Phi_File:%s Successed\n", fgm_info.point_phi_file_name);

				//用cuda做完傅里叶变换与逆变换后需要乘1/grid_numbers的系数来得到卷积的结果，因此将这步提前到初始化点电荷电势分布这里
				Scale_List << <ceilf((float)grid_numbers / 32), 32 >> >(grid_numbers, point_phi, grid_numbers_inverse);
				cufftExecR2C(plan_r2c, point_phi, fft_of_point_phi);
				printf("	Initial Point_Phi FFT Successed\n");
			}
		}
		//用有限差分法求解点电荷的电势分布
		else if ((fgm_info.fgm_mode & 0x000000F0) == CALCULATE_POINT_PHI_DAT)
		{
			//有限元计算点电荷分布的时候采用与分子模拟做FFT时同样的格点数目
			if ((fgm_info.fgm_mode & 0x00000F00) == SAME_FINITE_DIMENSION)
			{
				fd_info.grid_level = 1;
				fd_info.fd_grid_numbers = grid_numbers;
				fd_info.fd_grid_dimension = grid_dimension; 
				fd_info.fd_grid_length = grid_length;
				fd_info.fd_grid_length_inverse_square.x = 1. / grid_length.x / grid_length.x;
				fd_info.fd_grid_length_inverse_square.y = 1. / grid_length.y / grid_length.y;
				fd_info.fd_grid_length_inverse_square.z = 1. / grid_length.z / grid_length.z;
				fd_info.fd_neighbor = neighbor;
				printf("	Start Calculate Point_Phi\n");
				fd_info.Solving_Poisson_Equation(h_point_phi);
				printf("	Calculate Point_Phi Successed\n");
			}
			else
			{
				fd_info.fd_grid_numbers = fd_info.grid_level* fd_info.grid_level*fd_info.grid_level*grid_numbers;
				fd_info.fd_grid_dimension.int_x = fd_info.grid_level*grid_dimension.int_x;
				fd_info.fd_grid_dimension.int_y = fd_info.grid_level*grid_dimension.int_y;
				fd_info.fd_grid_dimension.int_z = fd_info.grid_level*grid_dimension.int_z;
				fd_info.fd_layer_numbers = fd_info.fd_grid_dimension.int_y*fd_info.fd_grid_dimension.int_x;
				fd_info.fd_grid_length.x = (float)grid_length.x / fd_info.grid_level;
				fd_info.fd_grid_length.y = (float)grid_length.y / fd_info.grid_level;
				fd_info.fd_grid_length.z = (float)grid_length.z / fd_info.grid_level;

				fd_info.fd_grid_length_inverse_square.x = 1. / fd_info.fd_grid_length.x / fd_info.fd_grid_length.x;
				fd_info.fd_grid_length_inverse_square.y = 1. / fd_info.fd_grid_length.y / fd_info.fd_grid_length.y;
				fd_info.fd_grid_length_inverse_square.z = 1. / fd_info.fd_grid_length.z / fd_info.fd_grid_length.z;
				fd_info.Initial_Grid_Neighbor();
				printf("	Start Calculate Point_Phi Level:%d\n", fd_info.grid_level);
				fd_info.Solving_Poisson_Equation(h_point_phi);
				printf("	Calculate Point_Phi Successed\n");
				free(fd_info.h_fd_neighbor), cudaFree(fd_info.fd_neighbor);

				//用cuda做完傅里叶变换与逆变换后需要乘1/grid_numbers的系数来得到卷积的结果，因此将这步提前到初始化点电荷电势分布这里
				cudaMemcpy(point_phi, h_point_phi, sizeof(float)*grid_numbers, cudaMemcpyHostToDevice);
				Scale_List << <ceilf((float)grid_numbers / 32), 32 >> >
					(grid_numbers, point_phi, grid_numbers_inverse);
				cufftExecR2C(plan_r2c, point_phi, fft_of_point_phi);
			}
		}
	}

	//格林函数积分的离散采样方法初始化
	if ((fgm_info.fgm_mode & 0x0000F000) == GREEN_SAMPLING_TEXTURE_VERSION)
	{
		Initial_TextureObj_Phi();
		gs_info.point_numbers = fgm_info.sphere_point_numbers;
		Malloc_Safely((void**)&gs_info.h_crd, sizeof(VECTOR)*gs_info.point_numbers);
		Cuda_Malloc_Safely((void**)&gs_info.crd, sizeof(VECTOR)*gs_info.point_numbers);
		Cuda_Malloc_Safely((void**)&gs_info.crd_with_ds, sizeof(VECTOR)*gs_info.point_numbers);
		Cuda_Malloc_Safely((void**)&gs_info.int_crd, sizeof(INT_VECTOR)*gs_info.point_numbers);

		FILE *sphere_point_crd_in = NULL;
		Open_File_Safely(&sphere_point_crd_in, fgm_info.sphere_pos_file_name, "rb");
		fread(gs_info.h_crd, sizeof(VECTOR), fgm_info.sphere_point_numbers, sphere_point_crd_in);
		fclose(sphere_point_crd_in);
		cudaMemcpy(gs_info.crd, gs_info.h_crd, sizeof(VECTOR)*gs_info.point_numbers, cudaMemcpyHostToDevice);
		cudaMemcpy(gs_info.crd_with_ds, gs_info.h_crd, sizeof(VECTOR)*gs_info.point_numbers, cudaMemcpyHostToDevice);

		gs_info.ds = (float)4.*const_parameter.pi * 3. / gs_info.sphere_radius / gs_info.sphere_radius / gs_info.point_numbers;
		Scale_List << <ceilf((float)3.*gs_info.point_numbers / 32), 32 >> >
			(3 * gs_info.point_numbers, (float*)gs_info.crd, gs_info.sphere_radius);
		Scale_List << <ceilf((float)3.*gs_info.point_numbers / 32), 32 >> >
			(3 * gs_info.point_numbers, (float*)gs_info.crd_with_ds, gs_info.sphere_radius*gs_info.ds);

		VECTOR crd_to_int_crd_cof = { 
			const_parameter.uint_max_float / fgm_info.box_length.x,
			const_parameter.uint_max_float / fgm_info.box_length.y,
			const_parameter.uint_max_float / fgm_info.box_length.z };
		Crd_To_Int_Crd << <ceilf((float)gs_info.point_numbers / 32), 32 >> >
			(gs_info.point_numbers, crd_to_int_crd_cof, gs_info.crd, gs_info.int_crd);

		printf("	Read Sphere_Point_File:%s Successed\n", fgm_info.point_phi_file_name);
	}
	printf("END (Initial Finite Green Method)\n");
	return 1;
}
int FINITE_GREEN_METHOD::Calculate_FGM_Force
(const int charge_numbers, const float *charge_density_list, 
const UNSIGNED_INT_VECTOR *uint_crd, VECTOR *frc)
{
	//清空电荷分布表
	Reset_List << <ceilf((float)grid_numbers / 128), 128 >> >(grid_numbers, rho, 0.);
	//用给的坐标与电荷密度计算得到电荷分布rho
	Distribute_Charge_To_Vertex << <ceilf((float)charge_numbers / 32), 32 >> >
		(charge_numbers, charge_density_list, uint_crd, uint_crd_to_grid_serial,
		grid_dimension, layer_numbers,
		neighbor, rho);
	//电荷分布rho做FFT与fft_of_point_phi相乘再逆FFT得到电势分布phi
	cufftExecR2C(plan_r2c, rho, fft_of_rho);
	cufft_times << <1, 512 >> >(grid_numbers, fft_of_point_phi, fft_of_rho, fft_of_phi);
	cufftExecC2R(plan_c2r, fft_of_phi, phi);//用cuda做完傅里叶变换与逆变换后需要乘1/grid_numbers的系数来得到卷积的结果，因此将这步提前到初始化点电荷电势分布这里

	//进行格林函数积分离散点采样
	if ((fgm_info.fgm_mode & 0x0000F000) == GREEN_SAMPLING_TEXTURE_VERSION)
	{
		cudaMemcpy3D(&copyParams_phi);
		Green_Force_Texture_Ver << <ceilf((float)charge_numbers / 128), 128 >> >
			(charge_numbers, uint_crd, charge_density_list,
			gs_info.point_numbers, gs_info.int_crd, gs_info.crd_with_ds,
			uint_crd_scale_to_1,
			half_grid_lenth_scale_to_1,
			texObj_phi, frc);
	}
	return 1;
}
int FINITE_GREEN_METHOD::Export_Phi(char *file_name)
{
	FILE *out = NULL;
	Open_File_Safely(&out, file_name, "wb");
	cudaMemcpy(h_phi, phi, sizeof(float)*grid_numbers, cudaMemcpyDeviceToHost);
	fwrite(h_phi, sizeof(float), grid_numbers, out);
	fclose(out);
	return 1;
}

//一些辅助函数
int Check_2357_Factor(int number)
{
	int tempn;
	while (number>0)
	{
		if (number == 1)
			return 1;
		tempn = number / 2;
		if (tempn * 2 != number)
			break;
		number = tempn;
	}

	while (number>0)
	{
		if (number == 1)
			return 1;
		tempn = number / 3;
		if (tempn * 3 != number)
			break;
		number = tempn;
	}

	while (number>0)
	{
		if (number == 1)
			return 1;
		tempn = number / 5;
		if (tempn * 5 != number)
			break;
		number = tempn;
	}

	while (number>0)
	{
		if (number == 1)
			return 1;
		tempn = number / 7;
		if (tempn * 7 != number)
			break;
		number = tempn;
	}

	return 0;
}
int Get_FFT_Patameter(float length, float factor)
{
	int tempi = (int)ceil(length / factor + 3) >> 2 << 2;

	while (1)
	{
		if (Check_2357_Factor(tempi))
			return tempi;
		tempi += 4;
	}
}
__global__ void cufft_times(int element_numbers, cufftComplex *a, cufftComplex *b, cufftComplex *c)
{
	for (int i = threadIdx.x; i < element_numbers; i = i + blockDim.x)
	{
		c[i] = cuCmulf(a[i], b[i]);
	}
}
__global__ void Distribute_Charge_To_Vertex
(const int atom_numbers, const float *charge_density, const UNSIGNED_INT_VECTOR *uint_crd,
const VECTOR scaler, const INT_VECTOR Nmax, const int Nxy,
const FINITE_GRID_NEIGHBOR *neighbor, float *rho
)
{
	int atom_i = blockDim.x*blockIdx.x + threadIdx.x;
	if (atom_i < atom_numbers)
	{

		float c_d = charge_density[atom_i];
		VECTOR lin;
		lin.x = (float)uint_crd[atom_i].uint_x*scaler.x;
		lin.y = (float)uint_crd[atom_i].uint_y*scaler.y;
		lin.z = (float)uint_crd[atom_i].uint_z*scaler.z;

		int grid_x = lin.x;
		int grid_y = lin.y;
		int grid_z = lin.z;

		grid_x = grid_x&((grid_x - Nmax.int_x) >> 31);//can opt
		grid_y = grid_y&((grid_y - Nmax.int_y) >> 31);
		grid_z = grid_z&((grid_z - Nmax.int_z) >> 31);
		int grid_serial = grid_z*Nxy + grid_y*Nmax.int_x + grid_x;

		lin.x = lin.x - (float)grid_x;
		lin.y = lin.y - (float)grid_y;
		lin.z = lin.z - (float)grid_z;

		float x = 1. - lin.x;
		float y = 1. - lin.y;
		float z = 1. - lin.z;

		float q0 = x*y*z*c_d;//need opt

		float qx = lin.x*y*z*c_d;
		float qy = x*lin.y*z*c_d;
		float qz = x*y*lin.z*c_d;
		float qxy = lin.x*lin.y*z*c_d;
		float qxz = lin.x*y*lin.z*c_d;
		float qyz = x*lin.y*lin.z*c_d;
		float qxyz = lin.x*lin.y*lin.z*c_d;

		atomicAdd(&rho[grid_serial], q0);
		int grid_serial_x = neighbor[grid_serial].x;
		int grid_serial_y = neighbor[grid_serial].y;
		int grid_serial_z = neighbor[grid_serial].z;
		atomicAdd(&rho[grid_serial_x], qx);
		atomicAdd(&rho[grid_serial_y], qy);
		atomicAdd(&rho[grid_serial_z], qz);

		grid_serial_z = neighbor[grid_serial_x].y;//need opt
		atomicAdd(&rho[grid_serial_z], qxy);
		grid_serial_z = neighbor[grid_serial_x].z;
		atomicAdd(&rho[grid_serial_z], qxz);
		grid_serial_z = neighbor[grid_serial_y].z;
		atomicAdd(&rho[grid_serial_z], qyz);
		grid_serial_z = neighbor[grid_serial_z].x;
		atomicAdd(&rho[grid_serial_z], qxyz);
	}
}
int FINITE_GREEN_METHOD::Initial_Grid_Neighbor()
{
	Cuda_Malloc_Safely((void**)&neighbor, sizeof(FINITE_GRID_NEIGHBOR)*grid_numbers);
	Malloc_Safely((void**)&h_neighbor, sizeof(FINITE_GRID_NEIGHBOR)*grid_numbers);
	int grid_serial;
	int grid_serial_2;
	int linx, liny, linz;
	for (int i = 0; i < grid_dimension.int_z; i = i + 1)
	{
		for (int j = 0; j < grid_dimension.int_y; j = j + 1)
		{
			for (int k = 0; k < grid_dimension.int_x; k = k + 1)
			{
				grid_serial = i*layer_numbers + j*grid_dimension.int_x + k;
				linx = k + 1;
				if (linx >= grid_dimension.int_x)
				{
					linx = linx - grid_dimension.int_x;
				}
				h_neighbor[grid_serial].x = i*layer_numbers + j*grid_dimension.int_x + linx;

				linx = k - 1;
				if (linx < 0)
				{
					linx = linx + grid_dimension.int_x;
				}
				h_neighbor[grid_serial].x_ = i*layer_numbers + j*grid_dimension.int_x + linx;

				liny = j + 1;
				if (liny >= grid_dimension.int_y)
				{
					liny = liny - grid_dimension.int_y;
				}
				h_neighbor[grid_serial].y = i*layer_numbers + liny*grid_dimension.int_x + k;

				liny = j - 1;
				if (liny < 0)
				{
					liny = liny + grid_dimension.int_y;
				}
				h_neighbor[grid_serial].y_ = i*layer_numbers + liny*grid_dimension.int_x + k;

				linz = i + 1;
				if (linz >= grid_dimension.int_z)
				{
					linz = linz - grid_dimension.int_z;
				}
				h_neighbor[grid_serial].z = linz*layer_numbers + j*grid_dimension.int_x + k;

				linz = i - 1;
				if (linz < 0)
				{
					linz = linz + grid_dimension.int_z;
				}
				h_neighbor[grid_serial].z_ = linz*layer_numbers + j*grid_dimension.int_x + k;

			}
		}
	}
	cudaMemcpy(neighbor, h_neighbor, sizeof(FINITE_GRID_NEIGHBOR)*grid_numbers, cudaMemcpyHostToDevice);
	return 1;
}
int FINITE_GREEN_METHOD::FINITE_DIFFERENCE_SOLVING_POISSON_EQUATION_INFORMATION::Solving_Poisson_Equation(float *h_point_phi)
{
	//进行一次泊松方程求解，考虑到格子数不一定与分子模拟过程中做fft的相同，因此求解中用的变量在这里声明并在这个函数内销毁
	float *rho = NULL;//体系的电荷分布
	float *h_rho = NULL;//体系的电荷分布
	float *h_point_phi_1 = NULL;//点电荷生成的电势分布
	float *point_phi_1 = NULL;//点电荷生成的电势分布
	float *point_phi_2 = NULL;//点电荷生成的电势分布
	float *point_phi_ref = NULL;//点电荷生成的电势分布用于收敛性判断
	Cuda_Malloc_Safely((void**)&rho, sizeof(float)* fd_grid_numbers);
	Malloc_Safely((void**)&h_rho, sizeof(float)* fd_grid_numbers);
	Malloc_Safely((void**)&h_point_phi_1, sizeof(float)* fd_grid_numbers);
	Cuda_Malloc_Safely((void**)&point_phi_1, sizeof(float)* fd_grid_numbers);
	Cuda_Malloc_Safely((void**)&point_phi_2, sizeof(float)* fd_grid_numbers);
	Cuda_Malloc_Safely((void**)&point_phi_ref, sizeof(float)* fd_grid_numbers);
	for (int i = 0; i < fd_grid_numbers; i = i + 1)
	{
		h_rho[i] = 0.;
	}
	h_rho[0] = 1. / fd_grid_length.x / fd_grid_length.y / fd_grid_length.z;//点电荷密度
	cudaMemcpy(rho, h_rho, sizeof(float)*fd_grid_numbers, cudaMemcpyHostToDevice);
	Reset_List << <ceilf((float)fd_grid_numbers / 128), 128 >> >(fd_grid_numbers, point_phi_1, 0.);
	Reset_List << <ceilf((float)fd_grid_numbers / 128), 128 >> >(fd_grid_numbers, point_phi_2, 0.);
	Reset_List << <ceilf((float)fd_grid_numbers / 128), 128 >> >(fd_grid_numbers, point_phi_ref, 1.);//参照phi一开始肯定应该不同于初始猜测的phi

	float h_phi_distance;
	float *phi_distance;//两次迭代的phi距离，表明收敛性
	Cuda_Malloc_Safely((void**)&phi_distance, sizeof(float));

	printf("	SOR: SOR_Step Convergence\n");
	Refresh_Gamma(first_gamma);
	for (int SOR_step = 0; SOR_step < first_iteration_steps; SOR_step = SOR_step + 1)
	{
		SOR_Iteration2 << <ceilf((float)fd_grid_numbers / 32), 32 >> >
			(fd_grid_numbers, fd_neighbor,
			rho, fd_grid_length_inverse_square,
			SOR_gamma_1, SOR_gamma_2, point_phi_1, point_phi_2);
		SOR_Iteration2 << <ceilf((float)fd_grid_numbers / 32), 32 >> >
			(fd_grid_numbers, fd_neighbor,
			rho, fd_grid_length_inverse_square,
			SOR_gamma_1, SOR_gamma_2, point_phi_2, point_phi_1);
		if (SOR_step % 100 == 0)
		{

			Calculate_Phi_Convergence << <1, 256 >> >
				(fd_grid_numbers, point_phi_1, point_phi_ref, phi_distance);
			Copy_List << <ceilf((float)fd_grid_numbers / 32), 32 >> >
				(fd_grid_numbers, point_phi_1, point_phi_ref);
			cudaMemcpy(&h_phi_distance, phi_distance, sizeof(float), cudaMemcpyDeviceToHost);
			printf("	    %8d %8E\n",SOR_step,(float)h_phi_distance / fd_grid_numbers);
		}
	}
	Refresh_Gamma(second_gamma);
	for (int SOR_step = 0; SOR_step < second_iteration_steps; SOR_step = SOR_step + 1)
	{
		SOR_Iteration2 << <ceilf((float)fd_grid_numbers / 32), 32 >> >
			(fd_grid_numbers, fd_neighbor,
			rho, fd_grid_length_inverse_square,
			SOR_gamma_1, SOR_gamma_2, point_phi_1, point_phi_2);
		SOR_Iteration2 << <ceilf((float)fd_grid_numbers / 32), 32 >> >
			(fd_grid_numbers, fd_neighbor,
			rho, fd_grid_length_inverse_square,
			SOR_gamma_1, SOR_gamma_2, point_phi_2, point_phi_1);
		if (SOR_step % 100 == 0)
		{

			Calculate_Phi_Convergence << <1, 256 >> >
				(fd_grid_numbers, point_phi_1, point_phi_ref, phi_distance);
			Copy_List << <ceilf((float)fd_grid_numbers / 32), 32 >> >
				(fd_grid_numbers, point_phi_1, point_phi_ref);
			cudaMemcpy(&h_phi_distance, phi_distance, sizeof(float), cudaMemcpyDeviceToHost);
			printf("	    %8d %8E\n", SOR_step, (float)h_phi_distance / fd_grid_numbers);
		}
	}

	SOR_Iteration2_Mean << <ceilf((float)fd_grid_numbers / 32), 32 >> >
		(fd_grid_numbers, point_phi_1, point_phi_2);
	cudaMemcpy(h_point_phi_1, point_phi_1, sizeof(float)*fd_grid_numbers, cudaMemcpyDeviceToHost);

	INT_VECTOR grid_dimension;
	grid_dimension.int_x = fd_grid_dimension.int_x / grid_level;
	grid_dimension.int_y = fd_grid_dimension.int_y / grid_level;
	grid_dimension.int_z = fd_grid_dimension.int_z / grid_level;
	int layer_numbers = grid_dimension.int_x*grid_dimension.int_y;
	for (int k = 0; k < grid_dimension.int_z; k = k + 1)
	{
		for (int j = 0; j < grid_dimension.int_y; j = j + 1)
		{
			for (int i = 0; i < grid_dimension.int_x; i = i + 1)
			{
				h_point_phi[k*layer_numbers + j*grid_dimension.int_x + i] = h_point_phi_1[grid_level*grid_level*grid_level*k*layer_numbers + grid_level*grid_level*j*grid_dimension.int_x + grid_level*i];
			}
		}
	}

	cudaFree(rho), cudaFree(point_phi_1), cudaFree(point_phi_2), cudaFree(point_phi_ref), cudaFree(phi_distance);
	free(h_rho), free(h_point_phi_1);
	return 1;
}
int FINITE_GREEN_METHOD::Set_Finite_Difference_Parameter(int grid_level, float first_gamma, int first_iteration_steps, float second_gamma, int second_iteration_steps)
{
	this->fd_info.grid_level = grid_level;
	this->fd_info.first_gamma = first_gamma;
	this->fd_info.first_iteration_steps = first_iteration_steps;
	this->fd_info.second_gamma = second_gamma;
	this->fd_info.second_iteration_steps = second_iteration_steps;
	return 1;
}
int FINITE_GREEN_METHOD::FINITE_DIFFERENCE_SOLVING_POISSON_EQUATION_INFORMATION::Refresh_Gamma(float target_gamma)
{
	SOR_gamma = target_gamma;
	SOR_gamma_1 = 1. - SOR_gamma;
	SOR_gamma_2 = 0.5*SOR_gamma /
		(fd_grid_length_inverse_square.x
		+ fd_grid_length_inverse_square.y
		+ fd_grid_length_inverse_square.z);
	return 1;
}
int FINITE_GREEN_METHOD::FINITE_DIFFERENCE_SOLVING_POISSON_EQUATION_INFORMATION::Initial_Grid_Neighbor()
{
	Cuda_Malloc_Safely((void**)&fd_neighbor, sizeof(FINITE_GRID_NEIGHBOR)*fd_grid_numbers);
	Malloc_Safely((void**)&h_fd_neighbor, sizeof(FINITE_GRID_NEIGHBOR)*fd_grid_numbers);
	int grid_serial;
	int grid_serial_2;
	int linx, liny, linz;
	for (int i = 0; i < fd_grid_dimension.int_z; i = i + 1)
	{
		for (int j = 0; j < fd_grid_dimension.int_y; j = j + 1)
		{
			for (int k = 0; k < fd_grid_dimension.int_x; k = k + 1)
			{
				grid_serial = i*fd_layer_numbers + j*fd_grid_dimension.int_x + k;
				linx = k + 1;
				if (linx >= fd_grid_dimension.int_x)
				{
					linx = linx - fd_grid_dimension.int_x;
				}
				h_fd_neighbor[grid_serial].x = i*fd_layer_numbers + j*fd_grid_dimension.int_x + linx;

				linx = k - 1;
				if (linx < 0)
				{
					linx = linx + fd_grid_dimension.int_x;
				}
				h_fd_neighbor[grid_serial].x_ = i*fd_layer_numbers + j*fd_grid_dimension.int_x + linx;

				liny = j + 1;
				if (liny >= fd_grid_dimension.int_y)
				{
					liny = liny - fd_grid_dimension.int_y;
				}
				h_fd_neighbor[grid_serial].y = i*fd_layer_numbers + liny*fd_grid_dimension.int_x + k;

				liny = j - 1;
				if (liny < 0)
				{
					liny = liny + fd_grid_dimension.int_y;
				}
				h_fd_neighbor[grid_serial].y_ = i*fd_layer_numbers + liny*fd_grid_dimension.int_x + k;

				linz = i + 1;
				if (linz >= fd_grid_dimension.int_z)
				{
					linz = linz - fd_grid_dimension.int_z;
				}
				h_fd_neighbor[grid_serial].z = linz*fd_layer_numbers + j*fd_grid_dimension.int_x + k;

				linz = i - 1;
				if (linz < 0)
				{
					linz = linz + fd_grid_dimension.int_z;
				}
				h_fd_neighbor[grid_serial].z_ = linz*fd_layer_numbers + j*fd_grid_dimension.int_x + k;

			}
		}
	}
	cudaMemcpy(fd_neighbor, h_fd_neighbor, sizeof(FINITE_GRID_NEIGHBOR)*fd_grid_numbers, cudaMemcpyHostToDevice);
	return 1;
}