#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define X 300   //X dimession of the data 
#define Y 2		//Y dimesnion of the data //TODO NEED TO KNOW THE DATA SIZE BEFORE IMPORTING
#define K 3		//NUMBER OF CLUSTERS TO DIVIDE THE DATA INTO
#define MAX_ITERS 10  //NUMBER OF ITERATIONS

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
		printf("COunts is %d \n", count);
		for (l=0;l<Y;l++){

			*(centroids+i*Y+l)=sum[l]/count;					
		}
	} 

}
void main(){

	FILE *fp;
	double num[X][Y];
	int    idx[X];
	double centroids[K][Y]={ {3,3}, {6,2}, {8,5} }; //initialization TODO make it random
	double num1;
	int i, j;
	fp=fopen("data.txt","r");
	if(fp==NULL) {
		printf("Exiting no file with such name \n");
		exit(-1);
	}
	printf("Hello World \n");
	for (i=0;i<X;i++){
		for (j=0;j<Y;j++){
			fscanf(fp,"%lf", &num1);
			num[i][j]=num1;
			//printf(" %.15lf ", num[i][j]);
		}
		//printf("\n");
	}
	fclose(fp);
	//printf("Double deref %e Single deref %p NO deref %p\n",**centroids, *centroids, centroids);

	for (i=0; i<MAX_ITERS; i++){

		findclosestcentroids((double *)num, (double *)centroids, &idx[0]);
			//for(int i=0; i<X; i++)  printf(" %d %d\n", i+1, idx[i]);  //Just for debugging
		computeCentroids((double *)num, &idx[0], (double *)centroids);

	}
	
	//Printing the final centroids
	for(i=0; i<K;i++){
			for(j=0; j<Y;j++){

				printf("%lf  ",centroids[i][j]);

			}
		printf("\n");
	}
	
}
