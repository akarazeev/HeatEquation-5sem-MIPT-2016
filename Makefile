CPP = /usr/local/opt/llvm/bin/clang
CPPFLAGS = -I/usr/local/opt/llvm/include -fopenmp
LDFLAGS = -L/usr/local/opt/llvm/lib

all:
	$(CPP) $(CPPFLAGS) main.c $(LDFLAGS)

#%.o: %.c %.h
#    $(CC) $(CFLAGS) -c $< -o $@
#
#CC = gcc
#CFLAGS = -Wall -Wextra --std=c99 -O3