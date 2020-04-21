#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<cuda.h>
#include<cuda_runtime.h>
#define X 300  //X dimession of the data 
#define Y 2		//Y dimesnion of the data //TODO NEED TO KNOW THE DATA SIZE BEFORE IMPORTING
#define K 3		//NUMBER OF CLUSTERS TO DIVIDE THE DATA INTO
#define MAX_ITERS 1  //NUMBER OF ITERATIONS

double *num=NULL;
double *centroids_c=NULL;
int *idx=NULL;

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
			printf("Distance of %d %e\n", index, sum);
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

__global__
void computeCentroids(double* num, int* idx, double* centroids_c){

	
	int index=blockIdx.x*blockDim.x+threadIdx.x;
	int stride=blockDim.x*gridDim.x;
	int offset=0; //offset keeps track if the same thread (number enters the loop the next time, as the thread id will be same)
	
	int i, m, j, l, count;
	double sum[Y]; 					//for(i=0;i<Y;i++) sum[i]=0.0;//is it reqd ?
	for(int i=index; i<K; i+=stride){

		int x=index+offset*stride;
		count=0;
		for(m=0;m<Y;m++) sum[m]=0.0;

		for(j =0; j<X; j++){

			if(idx[j]==x){

					count++;
					for (l=0;l<Y;l++){

						sum[l]=sum[l]+ *(num+j*Y+l);
					
					}
			
			}

		}
		printf("COunts is %d \n", count);
		for (l=0;l<Y;l++){

			*(centroids_c+x*Y+l)=sum[l]/count;					
		}
	}

}

int main(){

	FILE *fp;

	//initialization TODO make it random
	double centroids[K][Y]={{3,3},{6,2},{8,5}};

	double num1;
	int i, j, n_blocks, no_of_threads;
	
	no_of_threads=32;
	if (no_of_threads==(X*Y))
        n_blocks = (X*Y)/no_of_threads; //calculation of number of blocks baseds on threadscounts and world size
    else 
        n_blocks = (X*Y)/no_of_threads+1;;

	//Initializing CUDA memory 
	cudaMallocManaged(&num, sizeof(double)*X*Y);
	cudaMallocManaged(&centroids_c, sizeof(double)*K*Y);
	cudaMallocManaged(&idx, sizeof(int)*X);

	//Opening file and loading data into CUDA memory.
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
			//printf(" %.15lf ", num1);  //Just for debugging
		}
		//printf("\n");
	}
	fclose(fp);

	//Loading the 2-dimensional centroid Initialization data into a 1D VECTOR for CUDA. 
	for (i=0;i<K;i++){
		for (j=0;j<Y;j++){
			centroids_c[i*Y+j]=centroids[i][j];
			//printf(" %.15lf ", centroids_c[i*Y+j]);  //Just for debugging
		}
		//printf("\n");
	}
	for (i=0;i<MAX_ITERS;i++){

		findclosestcentroids<<< n_blocks, no_of_threads>>>(num, &centroids_c[0], &idx[0]);
		cudaDeviceSynchronize();
		computeCentroids<<<1,32>>>(num, &idx[0], &centroids_c[0]);
		cudaDeviceSynchronize();

	}
	
	for (i=0;i<X;i++){

		printf("%d===%d \n",i+1, idx[i]+1);
	}

	for(i=0; i<K;i++){
			for(j=0; j<Y;j++){

				printf("Centroids using CUDA %lf  ",*(centroids_c+i*Y+j));

			}
		printf("\n");
	}
	return 0;

}
