CC=gcc
INCDIRS=-I/usr/include/m4ri/
LDFLAGS=-lm4ri
override CFLAGS += -Wall $(INCDIRS) --std=c99 -O3
EXECUTABLES=main

all: $(EXECUTABLES)

%.o : %.c %.h
	$(CC) $(CFLAGS) -c $<

main.o: main.c isd.h utils.h simd.h
isd.o: isd.c isd.h utils.h simd.h
main: utils.o main.o isd.o
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

clean:
	rm -f $(EXECUTABLES) *.o
