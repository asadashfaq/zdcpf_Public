# Produced by CVXGEN, 2012-03-18 11:07:47 -0700.
# CVXGEN is Copyright (C) 2006-2011 Jacob Mattingley, jem@cvxgen.com.
# The code in this file is Copyright (C) 2006-2011 Jacob Mattingley.
# CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial
# applications without prior written permission from Jacob Mattingley.

# Filename: Makefile.
# Description: Basic Makefile.

OPT = -Wall

# for valgrind.
#OPT=-g -O0

OPT=-Os

# libmath is needed for sqrt. can remove in final product.
LDLIBS = -lm
CFLAGS = $(OPT) $(INCLUDES) -fpic
LIBFLAGS = -dynamiclib -WI,-soname
OBJS = solver.o matrix_support.o ldl.o uflowmin.o util.o
CC = gcc

.PHONY: all
all: testsolver

testsolver: solver.o matrix_support.o ldl.o testsolver.o util.o

# Include util.o for random functions and easy matrix printing.
#testsolver: solver.o matrix_support.o ldl.o util.o testsolver.o

solver.o: solver.h
matrix_support.o: solver.h
ldl.o: solver.h
util.o: solver.h
testsolver.o: solver.h

.PHONY : clean
clean :
	-rm -f *.o testsolver *.so

.PHONY: lib
lib: libuflowmin.dylib

libuflowmin.dylib: $(OBJS)
	$(CC) $(CFLAGS) $(LIBFLAGS),libuflowmin.dylib.1 -o libuflowmin.dylib $(OBJS)
