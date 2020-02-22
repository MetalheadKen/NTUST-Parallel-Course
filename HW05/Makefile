all: aos.exe soa.exe


.PHONY: all clean test

.DEFAULT_GOAL := all

CC=g++ 
CFLAGS=-O2 -Wall -march=native

stopwatch.o: stopwatch.cpp
	$(CC) $(CFLAGS) -c $< -o $@

YoUtil.o: YoUtil.cpp
	$(CC) $(CFLAGS) -c $< -o $@

%.o: hw05.cpp 
	$(CC) $(CFLAGS) -D$(basename $@) -c $< -o $@

%.exe: %.o stopwatch.o YoUtil.o
	$(CC) $(CFLAGS) -D$(basename $@) main.cpp -o $@ $(basename $@).o stopwatch.o YoUtil.o -lOpenCL

test: aos.exe soa.exe
	@echo Doing tests
	./aos.exe set1.pts 36 18 100 set1.ply
	./aos.exe set2.pts 36 18 100 set2.ply
	./aos.exe set3.pts 36 18 100 set3.ply
	./aos.exe set4.pts 36 18 100 set4.ply
	./aos.exe set5.pts 36 18 100 set5.ply
	./aos.exe set6.pts 36 18 100 set6.ply

clean: 
	@echo Cleanup artifacts 
	@-rm *.exe *.o *.ply >& /dev/null
