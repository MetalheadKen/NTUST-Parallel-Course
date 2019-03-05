all: main.exe

CC=g++
CFLAGS=-g -Wall

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
	valgrind ./main.exe check.pts M1.m check/AABB1.txt check/C1.txt check/o1.pts
	valgrind ./main.exe check.pts M2.m check/AABB2.txt check/C2.txt check/o2.pts
	valgrind ./main.exe check.pts M3.m check/AABB3.txt check/C3.txt check/o3.pts

# This cleans up the project
clean: 
	rm *.exe *.o
