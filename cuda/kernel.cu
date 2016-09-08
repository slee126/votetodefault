#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/sequence.h>

#ifndef __CUDACC__ 
#define __CUDACC__
#endif

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cuda.h>
#include <device_functions.h>
#include <cuda_runtime_api.h>

#include <iostream>
#include <stdio.h>

#define ARRAY_SIZE_X 32
#define ARRAY_SIZE_Y 16
#define ARRAY_SIZE_IN_BYTES ((ARRAY_SIZE_X) * (ARRAY_SIZE_Y) * (sizeof(unsigned int)))


__global__ void addKernel(const float *vo_init, const float *b, float *c, float *maxStore, float *globalMax, float *kgrid, float *bgrid, float *shockArray)
{
	//const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;



	__shared__ int N[1];
	N[0] = 0;

	__shared__ float maxVal[112];

	
	
	int rowIndex, colIndex;

	float beta, delta, alpha, ks, rf, sigma, theta, loss;
	float b0, k0, r, w, r_d, w_d, shock;

	beta = .95;
	sigma = 2;
	rf = .017;
	alpha = .36;
	theta = .28;
	loss = .99;
	delta = .025;

	ks = powf(((1.0 / beta - 1 + delta)*powf(1.0 / 3.0, alpha - 1) / alpha), 1 / (alpha - 1));

	if (threadIdx.x == 0){
		rowIndex = blockIdx.x / 224;
		colIndex = blockIdx.x % 224;

		k0 = .5*ks + rowIndex * .018318351020470;
		b0 = -1.5 + colIndex * .01;

		shock = shockArray[0];
		r = shock*alpha*powf(k0 * 3.0, alpha - 1);
		r_d = loss*shock;

		w = shock*(1-alpha)*powf(k0 * 3.0, alpha);
		w_d = loss*shock;



		kgrid[blockIdx.x] = rowIndex;
		bgrid[blockIdx.x] = colIndex;
		globalMax[blockIdx.x] = -1000.0;
	}

	__syncthreads();

	while (N[0] < 224){
		c[threadIdx.x + N[0] * 224] = vo_init[threadIdx.x + N[0] * 224] + b[threadIdx.x + N[0] * 224];
		__syncthreads();

		//Reduction *********************************************************************************
		if (threadIdx.x < 112){
			maxVal[threadIdx.x] = fmaxf(vo_init[threadIdx.x + N[0] * 224], vo_init[threadIdx.x + 112 + N[0] * 224]);
		}
		__syncthreads();
		if (threadIdx.x < 56){
			maxVal[threadIdx.x] = fmaxf(maxVal[threadIdx.x + 56], maxVal[threadIdx.x]);
		}
		__syncthreads();
		if (threadIdx.x < 28){
			maxVal[threadIdx.x] = fmaxf(maxVal[threadIdx.x + 28], maxVal[threadIdx.x]);
		}
		__syncthreads();
		if (threadIdx.x < 14){
			maxVal[threadIdx.x] = fmaxf(maxVal[threadIdx.x + 14], maxVal[threadIdx.x]);
		}
		__syncthreads();
		if (threadIdx.x < 8){
			maxVal[threadIdx.x] = fmaxf(maxVal[threadIdx.x + 8], maxVal[threadIdx.x]);
		}

		__syncthreads();
		if (threadIdx.x < 4){
			maxVal[threadIdx.x] = fmaxf(maxVal[threadIdx.x + 4], maxVal[threadIdx.x]);
		}

		__syncthreads();
		if (threadIdx.x < 2){
			maxVal[threadIdx.x] = fmaxf(maxVal[threadIdx.x + 2], maxVal[threadIdx.x]);
		}

		__syncthreads();
		if (threadIdx.x == 0){
			maxVal[threadIdx.x] = fmaxf(maxVal[threadIdx.x + 1], maxVal[threadIdx.x]);
			maxStore[N[0]] = maxVal[0];

			globalMax[blockIdx.x] = fmaxf(maxVal[0], globalMax[blockIdx.x]);
			N[0]++;
		}
		__syncthreads();
		//************************************************************************************************

	}//while


	__syncthreads();
}

__global__ void test(float *gpuDefault_decision, float *gpuPrice_mat, float *gpuShock, float *gpuProb, float *gpuVo_final, float *gpuVo_temp, float *gpuVnd_final,
	float *gpuVnd_temp, float *gpuVd_final, float *gpuVd_temp, float *gpuB_opt_nd, float *gpuK_opt_nd,  float *fpuK_opt_d)
{
	int rowIndex, colIndex;

	float beta, delta, alpha, ks, rf, sigma, theta, loss;
	float b0, k0, r, w, r_d, w_d, shock;
	float bp, kp, price, c;

	float tempDot=0.0;

	beta = .95;
	sigma = 2;
	rf = .017;
	alpha = .36;
	theta = .28;
	loss = .99;
	delta = .025;

	ks = powf(((1.0 / beta - 1 + delta)*powf(1.0 / 3.0, alpha - 1) / alpha), 1 / (alpha - 1));

	__shared__ float prob_vec[21];
	__shared__ float default_decision_vec[21];


	__shared__ int N[1];
	N[0] = 0;

	rowIndex = blockIdx.x / 224;
	colIndex = blockIdx.x % 224;
	k0 = .75*ks + rowIndex * .01;
	b0 = -.8 + colIndex * .005;
	

	__syncthreads();
	// This will do all the iterations for all the shocks
	for (int i = 0; i < 21; i++){
		shock = gpuShock[i];
		r = shock*alpha*powf(k0 * 3.0, alpha - 1);
		r_d = loss*shock;
		w = shock*(1 - alpha)*powf(k0 * 3.0, alpha);
		w_d = loss*shock;

		if (threadIdx.x == 0){
			//pull out  default decision
			for (int j = 0; j < 21; j++){
				prob_vec[j] = gpuProb[i * 21 + j];
				default_decision_vec[j] = gpuDefault_decision[rowIndex * 224 + colIndex + 224 * 224 * i];
				tempDot += prob_vec[j] * default_decision_vec[j];
				if (j == 20){
					gpuPrice_mat[blockIdx.x + i * 224 * 224] = (1 - tempDot) / (1 + .017);
					tempDot = 0;
				}
			}
		}
		__syncthreads();

		//entering thread independent
		while (N[0] < 224){
			bp = -.8 + threadIdx.x*.005;
			kp = .75*ks + N[0] * .01;
			price = gpuPrice_mat[i * 224 * 224 + N[0] * 224 + threadIdx.x];
			c = (1 + r - delta)*k0 + w / 3.0 - kp + price*bp + b0;

			__syncthreads();
			if (threadIdx.x ==0)
				N[0] += 1;
			__syncthreads();
		}//while N

	}//shocks i

}//function


int main()
{

	clock_t t, t1;
	t = clock();

	const  int blockNum = 50176;
	const int gridSize = 224;
	const int shockSize = 21;

	thrust::host_vector<float> shock(shockSize); 
	thrust::host_vector<float> prob(shockSize*shockSize); 

	thrust::host_vector<float> default_decision(gridSize * gridSize*shockSize);
	thrust::fill(default_decision.begin(), default_decision.end(), 1.0);

	thrust::host_vector<float> price_mat(gridSize * gridSize*shockSize);
	thrust::fill(price_mat.begin(), price_mat.end(), 1.0/(1.0+.017));

	thrust::host_vector<float> vo_final(gridSize * gridSize*shockSize);
	thrust::host_vector<float> vo_temp(gridSize * gridSize*shockSize);
	thrust::host_vector<float> vnd_final(gridSize * gridSize*shockSize);
	thrust::host_vector<float> vnd_temp(gridSize * gridSize*shockSize);

	thrust::host_vector<float> b_opt_nd(gridSize * gridSize*shockSize);
	thrust::host_vector<float> k_opt_nd(gridSize * gridSize*shockSize);
	
	thrust::host_vector<float> vd_temp(gridSize *shockSize);
	thrust::host_vector<float> vd_final(gridSize *shockSize);

	thrust::host_vector<float> k_opt_d(gridSize *shockSize);

	FILE *ptr_file;
	ptr_file = fopen("C:/Users/seung/Dropbox/matlab/vote/vo_final_init224.txt", "r");
	if (!ptr_file)
		printf("cant open");
	for (int i = 0; i < gridSize * gridSize * shockSize; i++)
		fscanf(ptr_file, "%f", &vo_final[i]);

	FILE *shock_file;
	shock_file = fopen("C:/Users/seung/Dropbox/matlab/vote/shocks.txt", "r");
	if (!shock_file)
		printf("cant open");
	for (int i = 0; i < shockSize; i++)
		fscanf(shock_file, "%f", &shock[i]);

	FILE *prob_file;
	prob_file = fopen("C:/Users/seung/Dropbox/matlab/vote/prob.txt", "r");
	if (!prob_file)
		printf("cant open");
	for (int i = 0; i < shockSize * shockSize; i++)
		fscanf(prob_file, "%f", &prob[i]);






	thrust::device_vector<float> gpuShock = shock;
	thrust::device_vector<float> gpuProb = prob;

	thrust::device_vector<float> gpuDefault_decision = default_decision;
	thrust::device_vector<float> gpuPrice_mat = price_mat;

	thrust::device_vector<float> gpuVo_final = vo_final;
	thrust::device_vector<float> gpuVo_temp = vo_final;
	thrust::device_vector<float> gpuVnd_final = vo_final;
	thrust::device_vector<float> gpuVnd_temp = vo_final;
	thrust::device_vector<float> gpuVd_final(gridSize*shockSize);
	thrust::device_vector<float> gpuVd_temp(gridSize*shockSize);
	thrust::device_vector<float> gpuB_opt_nd(gridSize*gridSize*shockSize);
	thrust::device_vector<float> gpuK_opt_nd(gridSize*gridSize*shockSize);
	thrust::device_vector<float> gpuK_opt_d(gridSize*shockSize);
	

	float* pgpuShock = thrust::raw_pointer_cast(&gpuShock[0]);
	float* pgpuProb = thrust::raw_pointer_cast(&gpuProb[0]);

	float* pgpuDefault_decision = thrust::raw_pointer_cast(&gpuDefault_decision[0]);
	float* pgpuPrice_mat = thrust::raw_pointer_cast(&gpuPrice_mat[0]);

	float* pgpuVo_final = thrust::raw_pointer_cast(&gpuVo_final[0]);
	float* pgpuVo_temp = thrust::raw_pointer_cast(&gpuVo_temp[0]);
	float* pgpuVnd_final = thrust::raw_pointer_cast(&gpuVnd_final[0]);
	float* pgpuVnd_temp = thrust::raw_pointer_cast(&gpuVnd_temp[0]);
	float* pgpuVd_final = thrust::raw_pointer_cast(&gpuVd_final[0]);
	float* pgpuVd_temp = thrust::raw_pointer_cast(&gpuVd_temp[0]);
	float* pgpuB_opt_nd = thrust::raw_pointer_cast(&gpuB_opt_nd[0]);
	float* pgpuK_opt_nd = thrust::raw_pointer_cast(&gpuK_opt_nd[0]);
	float* pgpuK_opt_d = thrust::raw_pointer_cast(&gpuK_opt_d[0]);

	test << <blockNum, gridSize >> >(pgpuDefault_decision, pgpuPrice_mat, pgpuShock, pgpuProb, pgpuVo_final, pgpuVo_temp, pgpuVnd_final, pgpuVnd_temp, pgpuVd_final, pgpuVd_temp, pgpuB_opt_nd, pgpuK_opt_nd, pgpuK_opt_d);

	//float* pgpuIn2 = thrust::raw_pointer_cast(&gpuIn1[0]);
	//float* pgpuIn1 = thrust::raw_pointer_cast(&gpuIn2[0]);
	//float* pgpuOut = thrust::raw_pointer_cast(&gpuOut[0]);
	//float* pgpuMax = thrust::raw_pointer_cast(&gpuMax[0]);
	//float* pgpuglobalMax = thrust::raw_pointer_cast(&gpuglobalMax[0]);

	//float *pgpuKgrid = thrust::raw_pointer_cast(&gpuKgrid[0]);
	//float *pgpuBgrid = thrust::raw_pointer_cast(&gpuBgrid[0]);

	//float *pgpuShock = thrust::raw_pointer_cast(&gpuShock[0]);


	////for (int k = 0; k < 100; k++){
	//int k = 0;
	//	addKernel << <blockNum, gridSize >> >(pgpuIn1, pgpuIn2, pgpuOut, pgpuMax, pgpuglobalMax, pgpuKgrid, pgpuBgrid, pgpuShock);

	//	thrust::host_vector<float> out = gpuOut; // initialize individual elements
	//	thrust::host_vector<float> out2 = gpuMax; // initialize individual elements
	//	thrust::host_vector<float> out3 = gpuglobalMax; // initialize individual elements

	//	thrust::host_vector<float> kgrid = gpuKgrid; // initialize individual elements
	//	thrust::host_vector<float> bgrid = gpuBgrid; // initialize individual elements

	//	gpuIn1 = gpuOut;

	//	t1 = clock() - t;
	//	std::cout << "Took " << (float)t1 / CLOCKS_PER_SEC << " for iter " << k << std::endl;
	////}

	//FILE *ptr_file2;
	//ptr_file2 = fopen("C:/Users/seung/Dropbox/matlab/vote/Added.txt", "w");

	//for (int h = 0; h < gridSize; h++){
	//	for (int i = 0; i < gridSize; i++){
	//		if (i < gridSize-1)
	//			fprintf(ptr_file2, "%f ", out[i + h * gridSize]);
	//		else
	//			fprintf(ptr_file2, "%f ", out[i + h * gridSize]);

	//	}
	//	fprintf(ptr_file2, "\n");
	//}


	//for (int i = 0; i < blockNum; i++)
	//	printf("for blockNum %d global max 1 is %f kgrid is %f bgrid is %f\n", i, out3[i], kgrid[i], bgrid[i]);



	t = clock() - t;
	std::cout << "Took " << (float)t / CLOCKS_PER_SEC << std::endl;
	std::cout << "DONE" << std::endl;

	int temp11;
	std::cin >> temp11;
	return 0;
}
