# C = /usr/local/opt/llvm/bin/clang
C = mpicc
# CFLAGS = -I/usr/local/opt/llvm/include -fopenmp
CFLAGS = -openmp
LDFLAGS = -L/usr/local/opt/llvm/lib

all:
	$(C) $(CFLAGS) ./src/main.c $(LDFLAGS) -o a.out

#%.o: %.c %.h
#    $(CC) $(CFLAGS) -c $< -o $@
#
#CC = gcc
#CFLAGS = -Wall -Wextra --std=c99 -O3
