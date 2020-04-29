#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<cuda.h>
#include<cuda_runtime.h>
#define X 16384  //X dimession of the data 
#define Y 3		//Y dimesnion of the data //TODO NEED TO KNOW THE DATA SIZE BEFORE IMPORTING
#define K 2		//NUMBER OF CLUSTERS TO DIVIDE THE DATA INTO
#define MAX_ITERS 5  //NUMBER OF ITERATIONS

double *num=NULL;
double *centroids_c=NULL;
int *idx=NULL;

__global__
void findclosestcentroids(double* num, double* centroids_c, int* idx){

	int index=blockIdx.x*blockDim.x+threadIdx.x;
	int stride=blockDim.x*gridDim.x;
	int offset=0; //offset keeps track if the same thread (number enters the loop the next time, as the thread id will be same)
	for(int i=index; i<X; i+=stride){
		
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

__global__
void computeCentroids(double* num, int* idx, double* centroids_c){

	
	int index=blockIdx.x*blockDim.x+threadIdx.x;
	int stride=blockDim.x*gridDim.x;
	int offset=0; //offset keeps track if the same thread (number enters the loop the next time, as the thread id will be same)
	
	int  m, j, l, count;
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
		//printf("Counts is %d \n", count);
		for (l=0;l<Y;l++){

			*(centroids_c+x*Y+l)=sum[l]/count;					
		}
	}

}

int main(){

	FILE *fp, *fw;

	int lower =0;
	int upper =X-1;
	srand(time(0));

	double num1;
	int i, j, k, n_blocks, no_of_threads,rnd_num;
	
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
	fp=fopen("input.txt","r");
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

	for (i = 0; i < K; i++) {

			rnd_num = (rand()%(upper-lower + 1)) + lower;
			//printf("%d ", rnd_num);  
			for (j=0;j<Y;j++){ 
        		*(centroids_c+i*Y+j) = *(num+rnd_num*Y+j);
        		//printf("Centroids are %e",*(centroids_c+i*Y+j)); 
        	} 
        //printf("\n");
    }

	for (i=0;i<MAX_ITERS;i++){

		int cudaDeviceCount;
		cudaError_t cE1,cE2;
		findclosestcentroids<<< n_blocks, no_of_threads>>>(num, centroids_c, idx);
		cE1=cudaGetDeviceCount( &cudaDeviceCount);
		cE2=cudaDeviceSynchronize();
		//printf("The rwo error is %d %d \n",cE1,cE2);
		const char* x_err=cudaGetErrorString (cE2);
		//printf("%s \n",x_err); 

		computeCentroids<<<1, 32>>>(num, &idx[0], centroids_c);
		cudaDeviceSynchronize();

	}

	/*
	for(i=0; i<K;i++){
			for(j=0; j<Y;j++){

				printf("Centroids using CUDA %lf  ",*(centroids_c+i*Y+j));

			}
		printf("\n");
	}
	*/

	for (i=0; i<X;i++){
		//printf("%d==%d\n",i+1, idx[i]+1);

		for (k=0;k<K;k++){

			if (idx[i]==k){

					for (j=0;j<Y;j++){			
						*(num+i*Y+j)=*(centroids_c+k*Y+j);
					}
			}
				
		}

	}
	fw=fopen("output.txt","w");
	
	for(i=0; i<X;i++){
	
		for(j=0; j<Y;j++){

				fprintf(fw,"%lf  ",*(num+i*Y+j));
				//	printf("%lf  ",num[i][j]);
			}
		fprintf(fw, "\n");
		//printf("\n");
	}

	fclose(fw);
	cudaFree(num);
	cudaFree(centroids_c);
	cudaFree(idx);
	return 0;

}
