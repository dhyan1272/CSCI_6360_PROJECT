#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#define X 786432   //X dimession of the data 
#define Y 3			//Y dimesnion of the data //TODO NEED TO KNOW THE DATA SIZE BEFORE IMPORTING
#define K 5			//NUMBER OF CLUSTERS TO DIVIDE THE DATA INTO
#define MAX_ITERS 5  //NUMBER OF ITERATIONS
double *num;
double* centroids;

void findclosestcentroids(double *num, double *centroids, int* idx){

	int i, j, l, min_ind; 
	double sum, dist[K],min_dist;
	for (i=0;i<X;i++){
		for (j=0;j<K;j++){
			sum=0;
			for (l=0;l<Y;l++){

				sum=sum+(*(num+i*Y+l)-*(centroids+j*Y+l))*(*(num+i*Y+l)-*(centroids+j*Y+l));
			}
			//printf("Distance %e\n",sum);
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
	*(idx+i)=min_ind;
	}
}

void computeCentroids(double* num, int* idx, double* centroids){

	int i,j, l, m, count;
	double sum[Y];
	for(i=0;i<Y;i++) sum[i]=0.0;

	for (i=0;i<K;i++){

		count=0;
		for(m=0;m<Y;m++) sum[m]=0.0;

		for(j =0; j<X; j++){

			if(idx[j]==i){

					count++;
					for (l=0;l<Y;l++){

						sum[l]=sum[l]+ *(num+j*Y+l);
					
					}
			
			}

		}
		//printf("COunts is %d \n", count);
		for (l=0;l<Y;l++){

			*(centroids+i*Y+l)=sum[l]/count;					
		}
	} 

}
void main(){

	FILE *fp, *fw;
	double num1;
	int i, j,k,rnd_num;
	num=(double*)calloc(X*Y, sizeof(double));
	centroids=(double*)calloc(K*Y,sizeof(double));
	int    idx[X];
	int lower=1;
	int upper =X;
	srand(time(0));
	fp=fopen("input.txt","r");
	if(fp==NULL) {
		printf("Exiting no file with such name \n");
		exit(-1);
	}
	printf("Hello World \n");
	for (i=0;i<X;i++){
		for (j=0;j<Y;j++){
			fscanf(fp,"%lf", &num1);
			*(num+i*Y+j)=num1;
		}
	}
	fclose(fp);

	for (i = 0; i < K; i++) {

			rnd_num = (rand()%(upper-lower + 1)) + lower;
			printf("%d ", rnd_num);  
			for (j=0;j<Y;j++){ 
        		*(centroids+i*Y+j) = *(num+rnd_num*Y+j); 
        	} 
        printf("\n");
    }

	//printf("Double deref %e Single deref %p NO deref %p\n",**centroids, *centroids, centroids);

	for (i=0; i<MAX_ITERS; i++){

		findclosestcentroids((double *)num, (double *)centroids, &idx[0]);
		computeCentroids((double *)num, &idx[0], (double *)centroids);

	}

	for(i=0; i<K;i++){
			for(j=0; j<Y;j++){

				printf("Centroids not using CUDA %lf  ",*(centroids+i*Y+j));

			}
		printf("\n");
	}

	for (i=0; i<X;i++){
		//printf("%d==%d\n",i+1, idx[i]+1);

		for (k=0;k<K;k++){

			if (idx[i]==k){

					for (j=0;j<Y;j++){			
						*(num+i*Y+j)=*(centroids+k*Y+j);
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
	free (num);
	free (centroids);
	
}