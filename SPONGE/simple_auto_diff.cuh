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

#ifndef SIMPLE_AUTO_DIFF_CUH
#define SIMPLE_AUTO_DIFF_CUH
#include "common.cuh"


/*XYJ备注：SAD=simple auto diff，简单自动微分
实现原理：利用操作符重载，将f(x,y)的关系同时用链式法则链接到df(x,y)上。效率肯定会有影响，暂时未明具体会影响多少
使用方法：1. 确定该部分需要求偏微分的数量，假设有1个，则后面使用的类就为SADfloat<1>，2个则为SADfloat<2>
2. 将包含微分的变量和过程用上面确定的类声明变量，其中对于想求的变量初始化时需要两个参数：本身的值和第i个变量
3. 正常计算，那么最后结果中的dval[i]即为第i个变量的微分。
使用样例：（均在No_PNC/generalized_Born.cu中）
1. 求有效伯恩半径对距离的导数：不求导数的函数为Effective_Born_Radii_Factor_CUDA，求导数的函数为GB_accumulate_Force_Energy_CUDA
2. 求GB能量对距离和有效伯恩半径的导数：不求导数的函数为GB_inej_Energy_CUDA，求导数的函数为GB_inej_Force_Energy_CUDA
*/
template<int N>
struct SADfloat
{
	float val;
	float dval[N];
	__device__ __host__ SADfloat<N>()
	{
		this->val = 0;
		for (int i = 0; i < N; i++)
		{
			this->dval[i] = 0;
		}
	}
	__device__ __host__ SADfloat<N>(int f, int id = -1)
	{
		this->val = (float)f;
		for (int i = 0; i < N; i++)
		{
			if (i != id)
				this->dval[i] = 0;
			else
				this->dval[i] = 1;
		}
	}
	__device__ __host__ SADfloat<N>(float f, int id = -1)
	{
		this->val = f;
		for (int i = 0; i < N; i++)
		{
			if (i != id)
				this->dval[i] = 0;
			else
				this->dval[i] = 1;
		}
	}
	__device__ __host__ SADfloat<N>(const SADfloat<N>& f)
	{
		this->val = f.val;
		for (int i = 0; i < N; i++)
		{
			this->dval[i] = f.dval[i];
		}
	}
	friend __device__ __host__ SADfloat<N> operator+ (const SADfloat<N>& f1, const SADfloat<N>& f2)
	{
		SADfloat<N> f;
		f.val = f1.val + f2.val;
		for (int i = 0; i < N; i++)
		{
			f.dval[i] = f1.dval[i] + f2.dval[i];
		}
		return f;
	}
	friend __device__ __host__ SADfloat<N> operator- (const SADfloat<N>& f1, const SADfloat<N>& f2)
	{
		SADfloat<N> f;
		f.val = f1.val - f2.val;
		for (int i = 0; i < N; i++)
		{
			f.dval[i] = f1.dval[i] - f2.dval[i];
		}
		return f;
	}
	friend __device__ __host__ SADfloat<N> operator* (const SADfloat<N>& f1, const SADfloat<N>& f2)
	{
		SADfloat<N> f;
		f.val = f1.val * f2.val;
		for (int i = 0; i < N; i++)
		{
			f.dval[i] = f2.val * f1.dval[i] + f1.val * f2.dval[i];
		}
		return f;
	}
	friend __device__ __host__ SADfloat<N> operator/ (const SADfloat<N>& f1, const SADfloat<N>& f2)
	{
		SADfloat<N> f;
		f.val = f1.val / f2.val;
		for (int i = 0; i < N; i++)
		{
			f.dval[i] = f1.dval[i] * f2.val - f2.dval[i] * f1.val;
			f.dval[i] /= f2.val * f2.val;
		}
		return f;
	}
	friend __device__ __host__ SADfloat<N> logf (const SADfloat<N>& f)
	{
		SADfloat<N> fa;
		fa.val = logf(f.val);
		for (int i = 0; i < N; i++)
		{
			fa.dval[i] = f.dval[i] / f.val;
		}

		return fa;
	}
	friend __device__ __host__ SADfloat<N> sqrtf(const SADfloat<N>& f)
	{
		SADfloat<N> fa;
		fa.val = sqrtf(f.val);
		for (int i = 0; i < N; i++)
		{
			fa.dval[i] = 0.5 / fa.val * f.dval[i];
		}
		return fa;
	}
	friend __device__ __host__ SADfloat<N> expf(const SADfloat<N>& f)
	{
		SADfloat<N> fa;
		fa.val = expf(f.val);
		for (int i = 0; i < N; i++)
		{
			fa.dval[i] = fa.val * f.dval[i];
		}
		return fa;
	}
};

#endif