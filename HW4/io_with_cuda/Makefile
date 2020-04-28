all: parallelio.c cuda.cu
	mpicc -g parallelio.c -c -o parallelio.o
	nvcc -g -G -arch=sm_70 cuda.cu -c -o parallelio-cuda.o
	mpicc -g parallelio.o parallelio-cuda.o -o parallelio-cuda-exe -L/usr/local/cuda-10.1/lib64/ -lcudadevrt -lcudart -lstdc++
