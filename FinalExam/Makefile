all: main.exe main-OMP.exe main-OCL.exe

CC=g++
#CC=clang++
CFLAGS=-std=c++11 -Wall -O2 -march=native -fomit-frame-pointer
OMP=-fopenmp 
#CFLAGS=-std=c++11 -g -march=native 
# CC=icpc
# CFLAGS=-std=c++11 -O2 -Wall -xHost -fomit-frame-pointer
# OMP=-qopenmp

lPNG=-lpng

# for detecting memory leaks
#PREFIX=valgrind --leak-check=full

.PRECIOUS: %.o

OBJS=PNGio.o stopwatch.o 

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@ 

FE-OMP.o: FE-OMP.cpp 
	$(CC) $(CFLAGS) $(OMP) -c $< -o $@ 

main.exe: main.cpp FE.o $(OBJS)
	$(CC) $(CFLAGS) $< -o $@ FE.o $(OBJS) $(lPNG)

main-OMP.exe: main.cpp FE-OMP.o $(OBJS)
	$(CC) $(CFLAGS) $(OMP) -D FEOMP $< -o $@ FE-OMP.o $(OBJS) $(lPNG)

main-OCL.exe: main.cpp FE-OCL.o YoUtil.o $(OBJS)
	$(CC) $(CFLAGS) -D FEOCL $< -o $@ FE-OCL.o $(OBJS) $(lPNG) YoUtil.o -lOpenCL

test: main.exe inputs/01.png inputs/02.png inputs/03.png inputs/04.png inputs/05.png inputs/06.png
	@ $(PREFIX) ./main.exe inputs/01.png outputs/out-01.png 80 110 1 80 6.0
	@ $(PREFIX) ./main.exe inputs/02.png outputs/out-02.png 6 64 1 128 7.2
	@ $(PREFIX) ./main.exe inputs/03.png outputs/out-03.png 4 70 2 80 4.8
	@ $(PREFIX) ./main.exe inputs/04.png outputs/out-04.png 180 200 5 160 5
	@ $(PREFIX) ./main.exe inputs/05.png outputs/out-05.png 15 20 1 32 3.2
	@ $(PREFIX) ./main.exe inputs/06.png outputs/out-06.png 40 110 1 80 6.0

testOMP: main-OMP.exe inputs/01.png inputs/02.png inputs/03.png inputs/04.png inputs/05.png inputs/06.png
	@ $(PREFIX) ./main-OMP.exe inputs/01.png outputs/out-01.png 80 110 1 80 6.0
	@ $(PREFIX) ./main-OMP.exe inputs/02.png outputs/out-02.png 6 64 1 128 7.2
	@ $(PREFIX) ./main-OMP.exe inputs/03.png outputs/out-03.png 4 70 2 80 4.8
	@ $(PREFIX) ./main-OMP.exe inputs/04.png outputs/out-04.png 180 200 5 160 5
	@ $(PREFIX) ./main-OMP.exe inputs/05.png outputs/out-05.png 15 20 1 32 3.2
	@ $(PREFIX) ./main-OMP.exe inputs/06.png outputs/out-06.png 40 110 1 80 6.0

testOCL: main-OCL.exe inputs/01.png inputs/02.png inputs/03.png inputs/04.png inputs/05.png inputs/06.png
	@ $(PREFIX) ./main-OCL.exe inputs/01.png outputs/out-01.png 80 110 1 80 6.0
	@ $(PREFIX) ./main-OCL.exe inputs/02.png outputs/out-02.png 6 64 1 128 7.2
	@ $(PREFIX) ./main-OCL.exe inputs/03.png outputs/out-03.png 4 70 2 80 4.8
	@ $(PREFIX) ./main-OCL.exe inputs/04.png outputs/out-04.png 180 200 5 160 5
	@ $(PREFIX) ./main-OCL.exe inputs/05.png outputs/out-05.png 15 20 1 32 3.2
	@ $(PREFIX) ./main-OCL.exe inputs/06.png outputs/out-06.png 40 110 1 80 6.0

testALL: test testOMP testOCL
clean: 
	rm *.o *.exe outputs/*
