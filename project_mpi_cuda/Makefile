all: gol-main.c gol-with-cuda.cu
	mpicc -g gol-main.c -c -o gol-main.o
	nvcc -g -G -arch=sm_70 gol-with-cuda.cu -c -o gol-with-cuda.o
	mpicc -g gol-main.o gol-with-cuda.o -o gol_all_exe -L /usr/local/cuda-10.1/lib64/ -lcudadevrt -lcudart -lstdc++