all: aos.exe soa.exe p1.exe p2.exe p3.exe p4.exe p5.exe p6.exe


.PHONY: all clean test


CC=g++ 
CFLAGS=-O2 -Wall -march=native

stopwatch.o: stopwatch.cpp
	$(CC) $(CFLAGS) -c $< -o $@

%.o: hw03.cpp 
	$(CC) $(CFLAGS) -D$(basename $@) -c $< -o $@

%.exe: %.o stopwatch.o
	$(CC) $(CFLAGS) -D$(basename $@) main.cpp -o $@ $(basename $@).o stopwatch.o 

test: aos.exe soa.exe p1.exe p2.exe p3.exe p4.exe p5.exe p6.exe
	@echo Doing tests
	./aos.exe set1.pts 18 9 100 100 set1.ply
	./soa.exe set1.pts 18 9 100 100 set1.ply

clean: 
	@echo Cleanup artifacts 
	@-rm *.exe *.o *.ply >& /dev/null
