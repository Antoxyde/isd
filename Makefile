CC=gcc
INCDIRS=-I/usr/include/m4ri/
LDFLAGS=-lm4ri
CFLAGS=-g -Wall -Werror $(INCDIRS) --std=c99 -O3
EXECUTABLES=test

all: $(EXECUTABLES)

test.o: test.c utils.h
test: utils.o test.o
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

clean:
	rm -f $(EXECUTABLES) *.o

%.o : %.c %.h
	$(CC) $(CFLAGS) -c $<


