#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>

int main(int argc, char *argv[]) {

	int myrank = 0;
    int numranks = 0;

	MPI_File fh;
	MPI_Status status;

	if( argc != 2 )
    {
		printf("The code requires the argument for the size of each blocks e.g. ./io 1024 \n");
		exit(-1);
    }

    int filesize = atoi(argv[1]);

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numranks);

    int* buf = malloc(filesize);

    int block_count = filesize/sizeof(int);

    for (int i=0; i<block_count; i++) { buf[i] = 1; }

    MPI_File_open(MPI_COMM_WORLD, "datafile.txt", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
 
	for (int i=0; i < 64; i++) {
		int write_index = filesize*myrank + filesize*numranks*i; 
		MPI_File_write_at(fh, write_index, buf, block_count, MPI_INT, &status);
	}

	MPI_Barrier(MPI_COMM_WORLD);

    MPI_File_close(&fh);

	MPI_Finalize();

	free(buf);

    return 1;
}
