all: main.exe

CC=g++
CFLAGS=-g -Wall -Werror -O3 -funroll-loops -march=native
VALGRIND=valgrind
VFLAGS=--leak-check=full --show-leak-kinds=all --trace-children=yes -q

OBJS=stopwatch.o hw01.o

%.o: %.cpp 
	$(CC) $(CFLAGS) -c $< -o $@ 

main.exe: $(OBJS) main.cpp
	$(CC) $(CFLAGS) main.cpp $(OBJS) -o $@

# This test use the check.pts, which contains 8 points that forms a cube.
# M1 shifts every point by (1.1, 2.2, 3.3)
# M2 scales every point by (1, 2, 3) 
# M3 rotates around Z-axis by 90 degrees
# valgrind is an utility that checks for memory leaks & other memory access issues.
test: main.exe
	$(VALGRIND) $(VFLAGS) ./main.exe check.pts M1.m check/AABB1.txt check/C1.txt check/o1.pts
	$(VALGRIND) $(VFLAGS) ./main.exe check.pts M2.m check/AABB2.txt check/C2.txt check/o2.pts
	$(VALGRIND) $(VFLAGS) ./main.exe check.pts M3.m check/AABB3.txt check/C3.txt check/o3.pts

SHELL := /bin/bash
LOOPS=10
PROG=./main.exe Site.pts M1.m check/AABB1.txt check/C1.txt check/o1.pts

loop: main.exe
	for ((i=1; i <= ${LOOPS}; ++i)) do $(PROG); done

# This cleans up the project
clean: 
	rm *.exe *.o
