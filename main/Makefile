CC = gcc
RM = rm -rf

USR_DIR = /home/9yelin9/.local

CFLAGS = -g -O2 -Wall -mcmodel=medium -I../include -I$(USR_DIR)/include -fopenmp
LDFLAGS = -L../lib -L$(USR_DIR)/lib -fopenmp 
LINKS = -ltrihb -lz -lm -lhdf5

OBJS = mc.o
TARGETS = ../mc

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

../mc : mc.o
	$(CC) $(LDLIBS) $(LDFLAGS) -o $@ mc.o $(LINKS)
