DEBUG ?= 1
ifeq ($(DEBUG), 1)
	CFLAGS=-g -Wextra -Werror -Wall $(INCDIRS) --std=c99 -Wl,-rpath=. -Wno-unused-function -mbmi
else
	CFLAGS=-march=native -Ofast $(INCDIRS) --std=c99 -Wl,-rpath=. -mbmi
endif


CC=gcc
INCDIRS=-I.
LDFLAGS=
EXECUTABLES=main

all: $(EXECUTABLES)

%.o : %.c %.h
	$(CC) $(CFLAGS) -c $<

main.o: main.c prange.h utils.h libpopcnt.h xoshiro256starstar.h
prange.o: prange.c prange.h utils.h libpopcnt.h xoshiro256starstar.h
table.o: table.c table.h utils.h
stern.o: stern.h stern.c utils.h
main: utils.o main.o prange.o libm4ri-0.0.20200125.so stern.o
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

clean:
	rm -f $(EXECUTABLES) *.o

