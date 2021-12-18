//#include "SITS_common.cuh"
//
//__device__ float log_add_log(float a, float b)
//{
//	return fmaxf(a, b) + logf(1.0 + expf(-fabsf(a - b)));
//}
//
//__global__ void SITS_Record_Ene(float *ene_record, const float *pw_ene, const float *pp_ene, const float pe_a, const float pe_b)
//{
//	float temp = *pw_ene * 0.5 + *pp_ene;
//	temp = pe_a * temp + pe_b;
//
//	*ene_record = temp;
////	printf("DEBUG ene_record: %f\n", ene_record[0]);
//}
//
//
//
//__global__ void SITS_Update_gf(const int kn, float *gf, 
//	const float *ene_record, const float *log_nk, const float *beta_k)
//{
//	int i = blockDim.x * blockIdx.x + threadIdx.x;
//	if (i < kn)
//	{
//		gf[i] = -beta_k[i] * ene_record[0] + log_nk[i];
////		printf("DEBUG gf: %d %f\n", i, gf[i]);
//	}
//}
//
//__global__ void SITS_Update_gfsum(const int kn, float *gfsum, const float *gf)
//{
//	if (threadIdx.x == 0)
//	{
//		gfsum[0] = -FLT_MAX;
//	};
//	for (int i = 0; i < kn; i = i + 1)
//	{
//		gfsum[0] = log_add_log(gfsum[0], gf[i]);
////		printf("DEBUG gfsum: %d %f %f\n", i, gfsum[0], gf[i]);
//	}
//}
//
//__global__ void SITS_Update_log_pk(const int kn, float *log_pk, 
//	const float *gf, const float *gfsum, const int reset)
//{
//	int i = blockDim.x * blockIdx.x + threadIdx.x;
//	if (i < kn)
//	{
//		float gfi = gf[i];
//		log_pk[i] = ((float)reset) * gfi + ((float)(1 - reset)) * log_add_log(log_pk[i], gfi - gfsum[0]);
//		//printf("DEBUG log_pk: %d %f %f\n", i, log_pk[i], gfsum[0]);
//	}
//}
//
//
//__global__ void SITS_Update_log_mk_inverse(const int kn, 
//	float *log_weight, float *log_mk_inverse, float *log_norm_old, 
//	float *log_norm, const float *log_pk, const float *log_nk)
//{
//	int i = blockDim.x * blockIdx.x + threadIdx.x;
//	if (i < kn - 1)
//	{
//		log_weight[i] = (log_pk[i] + log_pk[i + 1]) * 0.5;
//		//printf("DEBUG log_weight: %d %f %f\n", i, log_pk[i], log_pk[i + 1]);
//		log_mk_inverse[i] = log_nk[i] - log_nk[i + 1];
//		log_norm_old[i] = log_norm[i];
//		log_norm[i] = log_add_log(log_norm[i], log_weight[i]);
//		log_mk_inverse[i] = log_add_log(log_mk_inverse[i] + log_norm_old[i] - log_norm[i], log_pk[i + 1] - log_pk[i] + log_mk_inverse[i] + log_weight[i] - log_norm[i]);
//		//printf("DEBUG log_norm: %d %f %f\n", i, log_norm[i], log_weight[i]);
//	}
//}
//
//__global__ void SITS_Update_log_nk_inverse(const int kn,
//	float *log_nk_inverse, 	const float *log_mk_inverse)
//{
//	for (int i = 0; i < kn - 1; i++)
//	{
//		log_nk_inverse[i + 1] = log_nk_inverse[i] + log_mk_inverse[i];
//	}
//}
//
//__global__ void SITS_Update_nk(const int kn,
//	float *log_nk, float *nk, const float *log_nk_inverse)
//{
//	int i = blockDim.x * blockIdx.x + threadIdx.x;
//	if (i < kn )
//	{
//		log_nk[i] = -log_nk_inverse[i];
//		nk[i] = exp(log_nk[i]);
//	}
//}