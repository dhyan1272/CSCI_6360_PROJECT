#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<cuda.h>
#include<cuda_runtime.h>
#include<time.h>
#define X 300  //X dimession of the data 
#define Y 2		//Y dimesnion of the data //TODO NEED TO KNOW THE DATA SIZE BEFORE IMPORTING
#define K 3		//NUMBER OF CLUSTERS TO DIVIDE THE DATA INTO
#define MAX_ITERS 1  //NUMBER OF ITERATIONS
#define blockSize 32

double *num=NULL;
double *centroids_c=NULL;
int *idx=NULL;
double *sum=NULL;
int *count=NULL;

__global__
void findclosestcentroids(double* num, double* centroids_c, int* idx){

	int index=blockIdx.x*blockDim.x+threadIdx.x;
	int stride=blockDim.x*gridDim.x;
	int offset=0; //offset keeps track if the same thread (number enters the loop the next time, as the thread id will be same)
	for(int i=index; i<X*Y; i+=stride){
		
		int x=index+offset*stride;
		int j, l, min_ind; 
		double sum, dist[K],min_dist;
		
		for (j=0;j<K;j++){
			sum=0;
			for (l=0;l<Y;l++){
					sum=sum+(*(num+x*Y+l)-*(centroids_c+j*Y+l))*(*(num+x*Y+l)-*(centroids_c+j*Y+l));
			}
			dist[j]=sqrt(sum);
		}
		min_dist=dist[0];
		min_ind=0;
		for (j=0; j<K; j++){
			if (dist[j]<min_dist) {
				min_dist=dist[j];
				min_ind=j;
			}
		}
		*(idx+x)=min_ind;
		offset++;
	}
}

__global__ void reduce(double *g_idata, double *g_odata, int *g_odata_count, const unsigned int m, const int *idx, const unsigned int cl)
{
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x*(blockSize*2) + tid;
  if (idx[i]==cl) {
    unsigned int gridSize = blockSize*2*gridDim.x;
    const unsigned int n = X*Y;
    __shared__ double sdata[n];
    __shared__ double sdata_count[n];

    sdata[tid] = 0.0;
    sdata_count[tid] = 0;

    while (i < n) {
      sdata[tid] += g_idata[(i*Y+m)] + g_idata[(i*Y+m)+blockSize];
      sdata_count[tid] += 2; 
      i += gridSize;
    }
    __syncthreads();

    if (blockSize >= 512) {if (tid < 256) { sdata[tid] += sdata[tid + 256]; sdata_count[tid] += sdata_count[tid + 256] ;} __syncthreads(); }
    if (blockSize >= 256) {if (tid < 128) { sdata[tid] += sdata[tid + 128]; sdata_count[tid] += sdata_count[tid + 128]; } __syncthreads(); }
    if (blockSize >= 128) {if (tid < 64) { sdata[tid] += sdata[tid +   64]; sdata_count[tid] += sdata_count[tid +   64]; } __syncthreads(); }
    
    if (tid < 32) {
        if (blockSize >=  64) sdata[tid] += sdata[tid + 32];
        if (blockSize >=  32) sdata[tid] += sdata[tid + 16];
        if (blockSize >=  16) sdata[tid] += sdata[tid +  8];
        if (blockSize >=  8) sdata[tid] += sdata[tid +  4];
        if (blockSize >=  4) sdata[tid] += sdata[tid +  2];
        if (blockSize >=  2) sdata[tid] += sdata[tid +  1];

        if (blockSize >=  64) sdata_count[tid] += sdata_count[tid + 32];
        if (blockSize >=  32) sdata_count[tid] += sdata_count[tid + 16];
        if (blockSize >=  16) sdata_count[tid] += sdata_count[tid +  8];
        if (blockSize >=  8) sdata_count[tid] += sdata_count[tid +  4];
        if (blockSize >=  4) sdata_count[tid] += sdata_count[tid +  2];
        if (blockSize >=  2) sdata_count[tid] += sdata_count[tid +  1];
    }

    if (tid == 0) g_odata[blockIdx.x] = sdata[0];
    if (tid == 0) g_odata_count[blockIdx.x] = sdata_count[0];
  }
}

void computeCentroids(double* num, int* idx, double* centroids, int n_blocks) {
  double sum_g = 0.0;
  int count_g = 0;
  cudaMallocManaged(&sum, sizeof(double)*n_blocks);
  cudaMallocManaged(&count, sizeof(int)*n_blocks);
  for (int i=0; i<K; i++) {
    for (int m=0; m<Y; m++) {
      for(int j=0; j <n_blocks; ++j) {
	sum[j] = 0.0;
	count[j] = 0;
      }
      reduce <<<n_blocks, blockSize>>> (num, sum, count, m, idx, i);
      cudaDeviceSynchronize();

      //do blocksum
      sum_g = 0.0; count_g = 0;
      for(int j=0; j <n_blocks; ++j) {
	sum_g += sum[j];
	count_g += count[j];
      }
      printf("m = %d, sum=%f, count=%d \n", m, sum_g, count_g);
      *(centroids+i*Y+m)=sum_g/count_g;
    }
  }
}

int main(){
	/*To generate random centroids
        srand(time(0)); 
	for(int i = 0; i<K; i++) {
	  for(int i = 0; i<Y; i++) centroids[i][j] = rand()%10;
	}
	*/

	double centroids[K][Y]={{3,3},{6,2},{8,5}};
	double num1;
	int i, j, n_blocks;
	n_blocks = (X*Y + blockSize - 1)/blockSize;
/*
	if (no_of_threads==(X*Y))
        n_blocks = (X*Y)/no_of_threads; //calculation of number of blocks baseds on threadscounts and world size
    else 
        n_blocks = (X*Y)/no_of_threads+1;;
*/
	//Initializing CUDA memory 
	cudaMallocManaged(&num, sizeof(double)*X*Y);
	cudaMallocManaged(&centroids_c, sizeof(double)*K*Y);
	cudaMallocManaged(&idx, sizeof(int)*X);

	//Opening file and loading data into CUDA memory.
	FILE *fp;
	fp=fopen("data.txt","r");
	if(fp==NULL) {
		printf("Exiting no file with such name \n");
		exit(-1);
	}
	//Loading the 2-dimensional data into a 1D VECTOR for CUDA. 
	for (i=0;i<X;i++){
		for (j=0;j<Y;j++){
			fscanf(fp,"%lf", &num1);
			num[i*Y+j]=num1;
		}
	}
	fclose(fp);

	//Loading the 2-dimensional centroid Initialization data into a 1D VECTOR for CUDA. 
	for (i=0;i<K;i++){
		for (j=0;j<Y;j++){
			centroids_c[i*Y+j]=centroids[i][j];
		}
	}

	for (i=0;i<MAX_ITERS;i++){
//question: why are we passing the global arrays as arguments??, the functions already have the info
		findclosestcentroids <<<n_blocks, blockSize>>> (num, &centroids_c[0], &idx[0]);
		cudaDeviceSynchronize();
		computeCentroids(num, &idx[0], &centroids_c[0], n_blocks);
	}
	
	//for (i=0;i<X;i++){
	//	printf("%d===%d \n",i+1, idx[i]+1);
	//}
	return 0;
}
