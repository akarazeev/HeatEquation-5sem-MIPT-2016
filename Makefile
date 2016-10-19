all:
	gcc main.c

#%.o: %.c %.h
#    $(CC) $(CFLAGS) -c $< -o $@
#
#CC = gcc
#CFLAGS = -Wall -Wextra --std=c99 -O3