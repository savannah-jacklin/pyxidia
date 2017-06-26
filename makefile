
##
CC = gcc
CFLAGS=-I. -Wall
#CFLAGS=-I.
#LIBS = -lm
SCIAPPSLIB = /Users/sjacklin/Microlensing/cfitsio/
LIBS = -lm -L$(SCIAPPSLIB) -lcfitsio
#

all: simu

simu: simu.o 
	$(CC) $(CFLAGS) -o simu simu.o  $(LIBS)

clean:
	rm -f *.o simu 
