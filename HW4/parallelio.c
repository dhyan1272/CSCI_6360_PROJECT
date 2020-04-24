#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<unistd.h>

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

int main(int argc, char *argv[]) {

	unsigned long long write_start = 0;
	unsigned long long write_finish = 0;

	unsigned long long read_start = 0;
	unsigned long long read_finish = 0;

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

    if (myrank == 0) {
    	write_start = getticks();
    }

    for (int i=0; i<block_count; i++) { buf[i] = 1; }

    MPI_File_open(MPI_COMM_WORLD, "datafile.txt", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
 
	for (int i=0; i < 64; i++) {
		int write_index = filesize*myrank + filesize*numranks*i; 
		MPI_File_write_at(fh, write_index, buf, block_count, MPI_INT, &status);
	}

	MPI_Barrier(MPI_COMM_WORLD);

    MPI_File_close(&fh);

    if (myrank == 0) {
    	write_finish = getticks();
    	printf("Time taken to perform write operation: %llu ticks\n", (write_finish - write_start));
    }

    free(buf);

    buf = malloc(filesize);

    if (myrank == 0) {
    	read_start = getticks();
    }

  	MPI_File_open(MPI_COMM_WORLD, "datafile.txt", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  	MPI_File_set_view(fh, 0, MPI_INT, MPI_INT, "native", MPI_INFO_NULL);

	for (int i = 0; i < 64; i++) {
		int read_index = filesize*myrank + filesize*numranks*i; 
		MPI_File_read_at(fh, read_index, buf, block_count, MPI_INT, &status);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_File_close(&fh);

    if (myrank == 0) {
    	read_finish = getticks();
    	printf("Time taken to perform read operation: %llu ticks\n", (read_finish - read_start));
    }

	MPI_Finalize();

	free(buf);

    return 0;
}
