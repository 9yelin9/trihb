CC = gcc
RM = rm -rf

USR_DIR = /home/9yelin9/.local
#MPI_DIR = /opt/mpi/gcc-4.8.5/openmpi-4.1.0
MPI_DIR = /opt/mpi/gcc-7.3.1/openmpi-4.1.4/
HDF_DIR = /opt/hdf5/intel-2021.4.0/1.12.0

CFLAGS = -g -O2 -Wall -mcmodel=medium -I../include -I$(USR_DIR)/include -I$(MPI_DIR)/include -I$(HDF_DIR)/include -fopenmp
LDFLAGS = -L../lib -L$(USR_DIR)/lib -L$(MPI_DIR)/lib -L$(HDF_DIR)/lib -fopenmp 
LINKS = -ltrihb -lz -lm -lopenblas -lhdf5

OBJS = test.o mc.o
TARGETS = ../test ../mc

.PHONY : all clean dep
.SUFFIXES : .c .o

.c .o :
	$(CC) $(CFLAGS) -c $<

all : $(TARGETS)
clean :
	$(RM) *.o
	$(RM) $(TARGETS)
dep :
	$(CC) $(CFLAGS) -M $(OBJS:.o=.c) 

../test : test.o
	$(CC) $(LDLIBS) $(LDFLAGS) -o $@ test.o $(LINKS)
../mc : mc.o
	$(CC) $(LDLIBS) $(LDFLAGS) -o $@ mc.o $(LINKS)
../obs : obs.o
	$(CC) $(LDLIBS) $(LDFLAGS) -o $@ obs.o $(LINKS)
