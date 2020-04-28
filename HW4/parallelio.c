/*
 * The code has been written by Karan Bhanot, Abolaji Adesoji, Aditya Joshi and Dhyanjyoti Nath.
 *
 * Some function definitions are referenced from
 * sample code provided by Christopher D. Carothers, 
 * provided as part of his class assignment of Parallel Computing 
 * Spring 2020.
 */

// Import required libraries
#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<unistd.h>

// Creating an alias to track ticks
typedef unsigned long long ticks;

/*
 * Function that returns the ticks at the momment in time
 * (Taken from the code provided as part of Parallel Computing class)
 */
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

/*
 * main function that handles the MPI operations for writing and reading
 * a file using parallel IO while keeping track of the number of ticks
 */
int main(int argc, char *argv[]) {

  // Variables to track ticks for writing and reading
  unsigned long long write_start = 0.0;
  unsigned long long write_finish = 0.0;
  unsigned long long read_start = 0.0;
  unsigned long long read_finish = 0.0;

  // Variables to track MPI rank and total ranks
  int myrank = 0;
  int numranks = 0;
  int count = 0;

  // MPI variables for file and status of reading
  MPI_File fh;
  MPI_Status status;

  // Check if the required number of arguments were provided
  if( argc != 2 )
  {
    printf("The code requires the argument for the size of each blocks e.g. ./parallelio-cuda-exe 1024 \n");
    exit(-1);
  }

  // Assign value based on the provided value
  long long blocksize = atoll(argv[1]);

  // Initialize MPI
  MPI_Init(&argc, &argv);

  // Get current rank and total number of ranks
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &numranks);

  // Allocate memory to buffer
  long long* buf = malloc(blocksize);

  // Calculate the total count of values in each block to write/read
  long long block_count = blocksize/sizeof(long long);

  // Set each value to be 1
  for (long long i=0; i<block_count; i++) { buf[i] = 1; }

  /*************************************************
   **************************************************
   ***************** WRITING TO FILE ****************
   **************************************************
   *************************************************/

  // If rank is 0, start calculation using ticks
  if (myrank == 0) {
    write_start = getticks();
  }

  // Open the file in write mode. If it does not exist, create it first
  MPI_File_open(MPI_COMM_WORLD, "datafile.txt", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

  // Write the data for each block in accordance with the rank number
  // If there are 64 ranks, the order of writing is:
  // Rank 0 - Block 0
  // Rank 1 - Block 0
  // .
  // .
  // .
  // Rank 64 - Block 0
  // Rank 0 - Block 1
  // .
  // .
  for (int i=0; i < 64; i++) {
    long long write_index = blocksize*myrank + blocksize*numranks*i; 
    MPI_File_write_at(fh, write_index, buf, block_count, MPI_LONG_LONG, &status);
    MPI_Get_count(&status, MPI_LONG_LONG, &count);
    if (count != block_count) {
        printf("Error writing data to file.\n");	
    }
}
  // Add barrier
  MPI_Barrier(MPI_COMM_WORLD);

  // After all blocks and rank writing, when rank is 0, print the total ticks
  if (myrank == 0) {
    write_finish = getticks();
    printf("Time taken to perform write operation: %.3f seconds\n", (write_finish - write_start)/52000000.0f);
  }

  // Close the file
  MPI_File_close(&fh);

  // Free the buffer
  free(buf);

  // Reinitialize buffer
  buf = malloc(blocksize);

  /*************************************************
   **************************************************
   *************** READING FROM FILE ****************
   **************************************************
   *************************************************/

  // If rank is 0, start calculation using ticks
  if (myrank == 0) {
    read_start = getticks();
  }

  // Open the file in write mode. If it does not exist, create it first
  MPI_File_open(MPI_COMM_WORLD, "datafile.txt", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

  // Read in all the data using various ranks and blocks
  // If there are 64 ranks, the order of reading is:
  // Rank 0 - Block 0
  // Rank 1 - Block 0
  // .
  // .
  // .
  // Rank 64 - Block 0
  // Rank 0 - Block 1
  // .
  // .
  for (int i = 0; i < 64; i++) {
    long long read_index = blocksize*myrank + blocksize*numranks*i; 
    MPI_File_read_at(fh, read_index, buf, block_count, MPI_LONG_LONG, &status);

    MPI_Get_count(&status, MPI_LONG_LONG, &count);
    if (count != block_count) {
        printf("Error reading data from file.\n");
    }
  }

  // Add barrier
  MPI_Barrier(MPI_COMM_WORLD);

  // After all blocks and rank writing, when rank is 0, print the total ticks
  if (myrank == 0) {
    read_finish = getticks();
    printf("Time taken to perform read operation: %.3f seconds\n", (read_finish - read_start)/52000000.0f);
  }

  //Close the file
  MPI_File_close(&fh);

  // Free the buffer
  free(buf);

  // Finalize MPI
  MPI_Finalize();

  return 0;
}
