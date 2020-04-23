#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<stdbool.h>


int main(int argc, char *argv[])
{
  int rank = 0;
  int nprocs = 0;
  double start_time, end_time;

  printf("This is the Assignment 4.\n");

  if( argc != 2 )
  {
    printf("It requires 1 argument: Block Size \n");
    exit(-1);
  }

  // MPI Setup
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
  MPI_Status status;
  MPI_Request request1, request2, request3, request4;

  MPI_File fh;
  int count, errs;
  int num_block = 64;
  int block_size = atoi(argv[1]);
  int block_count = (block_size)/sizeof(int);
  int* buf = malloc(block_count *sizeof(int));
  for (int i = 0; i <block_count; i++) {
    buf[i] = 1;
  }

  MPI_File_open(MPI_COMM_WORLD, "~/scratch/filename.txt", MPI_MODE_CREATE | MPI_MODE_RDWR | MPI_MODE_DELETE_ON_CLOSE, MPI_INFO_NULL, &fh);
  MPI_File_set_view(fh, 0, MPI_INT, MPI_INT, "native", MPI_INFO_NULL);

  for (int iny = 0; iny < num_block; iny++) {
  
    int block_start = (block_count * rank) + (((iny + num_block)%num_block) * (nprocs*block_count));
    MPI_File_write_at(fh, block_start, buf, block_count, MPI_INT, &status);
    MPI_Get_count( &status, MPI_INT, &count );
    // confirm the number of ints written
    if (count != block_count) { //confirm this
        errs++;
        fprintf( stderr, "Wrong count (%d) on write-ordered\n", count );fflush(stderr);
    }
  }
  MPI_File_close(&fh);

  MPI_File_open(MPI_COMM_WORLD, "~/scratch/filename.txt", MPI_MODE_CREATE | MPI_MODE_RDWR | MPI_MODE_DELETE_ON_CLOSE, MPI_INFO_NULL, &fh);
  MPI_File_set_view(fh, 0, MPI_INT, MPI_INT, "native", MPI_INFO_NULL);

  for (int iny = 0; iny < num_block; iny++) {
    int block_start = (block_count * rank) + (((iny + num_block)%num_block) * (nprocs*block_count));
    MPI_File_read_at(fh, block_start, buf, block_count, MPI_INT, &status);
  }
  MPI_File_close(&fh);

  MPI_Finalize();
  free(buf);

  return true;
}
