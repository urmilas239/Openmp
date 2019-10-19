# Adapted from https://sites.google.com/lbl.gov/cs267-spr2018/hw-2?authuser=0

CC = g++
MPCC = mpicxx
OPENMP = -fopenmp 
CFLAGS = -O3 -std=c++11
LIBS = 


TARGETS = serial binned binnedopenmp openmp mpi autograder

all:	$(TARGETS)

serial: serial.o common.o
	$(CC) -o $@ $(LIBS) serial.o common.o -lm
binned: binned.o common.o
	$(CC) -o $@ $(LIBS) binned.o common.o -lm
binnedopenmp: binnedopenmp.o common.o
	$(CC) -o $@ $(LIBS) $(OPENMP)  binnedopenmp.o common.o -lm 
autograder: autograder.o common.o
	$(CC) -o $@ $(LIBS) autograder.o common.o -lm
openmp: openmp.o common.o
	$(CC) -o $@ $(LIBS) $(OPENMP) openmp.o common.o -lm 
mpi: mpi.o common.o
	$(MPCC) -o $@ $(LIBS) $(MPILIBS) mpi.o common.o -lm

autograder.o: autograder.cpp common.h
	$(CC) -c $(CFLAGS) autograder.cpp 
openmp.o: openmp.cpp common.h
	$(CC) -c $(OPENMP) $(CFLAGS) openmp.cpp
binnedopenmp.o: binnedopenmp.cpp common.h
	$(CC) -c $(OPENMP) $(CFLAGS) binnedopenmp.cpp
serial.o: serial.cpp common.h
	$(CC) -c $(CFLAGS) serial.cpp
binned.o: binned.cpp common.h
	$(CC) -c $(CFLAGS) binned.cpp
mpi.o: mpi.cpp common.h
	$(MPCC) -c $(CFLAGS) mpi.cpp
common.o: common.cpp common.h
	$(CC) -c  $(CFLAGS) common.cpp

clean:
	rm -f *.o $(TARGETS) *.stdout *.txt *.error
