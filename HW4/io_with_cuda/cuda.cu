/*
* The code has been written by Karan Bhanot, Abolaji Adesoji, Aditya Joshi and Dhyanjyoti Nath.
*
* Some function definitions are referenced from
* sample code provided by Christopher D. Carothers, 
* provided as part of his class assignment of Parallel Computing 
* Spring 2020.
*/

// Include headers (including CUDA)
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<stdbool.h>
#include<cuda.h>
#include<cuda_runtime.h>

// Buffer
extern long long *buf;

/*
* Returns the inialized buffer on CUDA
*/
extern "C" void getBuffer( int rank, int numranks, long long blocksize )
{
    // Check and assign the device for this MPI rank
	cudaError_t cE;
	int cudaDeviceCount;

    // Check if enough devices are available
	if ((cE = cudaGetDeviceCount(&cudaDeviceCount)) != cudaSuccess) {
		printf("Unable to determine cuda device count, error is %d, count is %d\n", cE, cudaDeviceCount);
        exit(-1);
	}

    // Set device given that it is available
	if ((cE = cudaSetDevice(rank % cudaDeviceCount)) != cudaSuccess) {
		printf(" Unable to have rank %d set to cuda device %d, error is %d \n", rank, (rank % cudaDeviceCount), cE);
        exit(-1);
	}
	// Assign memory to the buf variable
	cudaMallocManaged(&buf, blocksize);
}
