#include<stdio.h>
#include<iostream>
#include<string.h>
#include<mpi.h>
//#include <netcdf>

int main(int argc, char* argv[]){
  
  int myrank, numranks,y;
  std::string fname;
  //fname=argv[1];
  MPI_File fh;
  MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &numranks);
  char const* x="Hello";
  char xwrite[15];
  y= strlen(x);

  MPI_File_open(MPI_COMM_WORLD, "TEST1.txt", MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &fh);
  MPI_File_write_at(fh,myrank, &x[myrank%y], sizeof(char), MPI_CHAR, &status);
  MPI_File_read_at(fh, myrank, &xwrite[myrank], sizeof(char), MPI_CHAR, &status);
  printf("For rank %d, the value read is %c \n", myrank, xwrite[myrank]);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_File_close(&fh);
  MPI_Finalize();

}
