#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<stdbool.h>
#include<time.h>
#include<mpi.h>
#include<string.h>
#include<math.h>

#define X 1049088			//X dimession of the data 
#define Y 3				//Y dimesnion of the data //TODO NEED TO KNOW THE DATA SIZE BEFORE IMPORTING
#define K 3				//NUMBER OF CLUSTERS TO DIVIDE THE DATA INTO
#define MAX_ITERS 10  	//NUMBER OF ITERATIONS

double *num;
double *centroids_c;
double *centroids_cresult;
int *idx;

void k_means_kernel_launch(double*, double*, int*, int, int, int);
void cuda_init(int, int, int);
void cuda_free(double*, double*, double*, int*);
void assign(double*, double*, int*, int, int, int);

typedef unsigned long long ticks;

static __inline__ ticks getticks(void)
{
  unsigned int tbl, tbu0, tbu1;

  do {
    __asm__ __volatile__ ("mftbu %0" : "=r"(tbu0));
    __asm__ __volatile__ ("mftb %0" : "=r"(tbl));
    __asm__ __volatile__ ("mftbu %0" : "=r"(tbu1));
  } while (tbu0 != tbu1);

  return (((unsigned long long)tbu0) << 32) | tbl;
}

int main(int argc, char *argv[]){

	// Define variables to keep track of time
	unsigned long long start = 0;
	unsigned long long finish = 0;

	int myrank, numranks, result;
	int i,each_chunk,j,k, no_of_threads,n_blocks;
    double starttime, endtime;
   
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numranks);
    MPI_Request request1, request2, request3, request4; 
    MPI_Status status;
    MPI_File fh;

    each_chunk=X/numranks;
	if(myrank==numranks-1)
		each_chunk=each_chunk+X%numranks;

	no_of_threads=32;
	if (no_of_threads==each_chunk)
        n_blocks = each_chunk/no_of_threads; //calculation of number of blocks baseds on threadscounts and world size
    else 
        n_blocks = each_chunk/no_of_threads+1;

    cuda_init(each_chunk, myrank, numranks);

	result=MPI_File_open(MPI_COMM_WORLD, "input.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
		if(result != MPI_SUCCESS) {printf("Error in opening the file\n"); exit(-1);}

	result=MPI_File_read_at(fh, myrank*each_chunk*Y*sizeof(double), num, each_chunk*Y, MPI_DOUBLE, &status);
		if(result != MPI_SUCCESS) {printf("Error in reading the file\n"); exit(-1);}
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_File_close(&fh);

	// Starting K-Means clustering

	if (myrank==0){

		// Let us start calculation of time
		start = getticks();

		int lower =0;
		int upper =each_chunk-1;
		srand(time(0));
		for (int i = 0; i < K; i++) {

			int rnd_num = (rand()%(upper-lower + 1)) + lower;
			//printf("%d ", rnd_num);  
			
			for (int j=0;j<Y;j++){ 
        		*(centroids_c+i*Y+j) = *(num+rnd_num*Y+j);
        		//printf("Centroids are %e",*(centroids_c+i*Y+j)); 
        	} 
        //printf("\n");
    	}

	}

	MPI_Bcast(centroids_c, K*Y, MPI_DOUBLE, 0, MPI_COMM_WORLD);
 
    for (i=0;i<MAX_ITERS;i++){


   		k_means_kernel_launch(num, centroids_c,idx,each_chunk, n_blocks, no_of_threads);
	    MPI_Barrier(MPI_COMM_WORLD);
	    MPI_Allreduce(centroids_c, centroids_cresult, K*Y, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	   	for (int a=0; a<K;a++){
			
			for (int b=0;b<Y;b++){

					*(centroids_c+a*Y+b)=*(centroids_cresult+a*Y+b)/numranks;
	
				}
	
		}


	    //double* temp;
	    //temp=centroids_cresult;
	    //centroids_cresult=centroids_c;
	    //centroids_cresult=temp;
    }
    
	//if(!myrank)
    for (i=0; i<K;i++){
		
		for (k=0;k<Y;k++){

			printf ("Centroids %d new %e", myrank, *(centroids_c+i*Y+k));
		}
		printf("\n");
	}
    	
   
    /*

    for (i=0; i<each_chunk;i++){
		//printf("%d==%d\n",i+1, idx[i]+1);

		for (k=0;k<K;k++){

			if (idx[i]==k){

					for (j=0;j<Y;j++){			
						*(num+i*Y+j)=*(centroids_cresult+k*Y+j)/numranks;
					}
			}
				
		}

	}
	*/
	//Parallelizing the above operation.
	assign(num, centroids_c, idx, each_chunk, n_blocks, no_of_threads);

	MPI_Barrier(MPI_COMM_WORLD);

	// KMeans algorithm completed
	
	if (myrank == 0) {
		// Close timer and print total time taken
		finish = getticks();
		printf("Total time taken to run K-Means on %d pixel image with %d clusters is %llu seconds.\n", X, K, (finish-start)/512000000.0f);
	}


   	result=MPI_File_open(MPI_COMM_WORLD, "output.bin",  MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &fh);
		if(result != MPI_SUCCESS) {printf("Error in opening the file\n"); exit(-1);}

	result=MPI_File_write_at(fh, myrank*each_chunk*Y*sizeof(double), num, each_chunk*Y, MPI_DOUBLE, &status);
		if(result != MPI_SUCCESS) {printf("Error in writing the file\n"); exit(-1);}


	MPI_Barrier(MPI_COMM_WORLD);
	MPI_File_close(&fh);

	cuda_free(num, centroids_c, centroids_cresult, idx);

    MPI_Finalize();
	return 0;
}