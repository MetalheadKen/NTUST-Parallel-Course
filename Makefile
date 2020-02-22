all: main.exe

CC=g++
CFLAGS=-Wall -g -march=native -O3

%.exe: %.cpp
	$(CC) $(CFLAGS) $< -o $@


OBJS=stopwatch.o hw02.o

%.o: %.cpp %.hpp
	$(CC) $(CFLAGS) -c $< -o $@

main.exe: main.cpp $(OBJS)
	$(CC) $(CFLAGS) $< -o $@ $(OBJS)

test: main.exe set1.pts set2.pts set3.pts set4.pts set5.pts set6.pts
	valgrind -q ./$< set1.pts 18 9 100 100 set1.ply
	valgrind -q ./$< set2.pts 10 10 10 20 set2.ply
	valgrind -q ./$< set3.pts 10 10 10 20 set3.ply
	valgrind -q ./$< set4.pts 10 10 10 20 set4.ply
	valgrind -q ./$< set5.pts 10 10 10 20 set5.ply
	valgrind -q ./$< set6.pts 10 10 10 20 set6.ply

clean: 
	rm set*.ply *.exe *.o
	
