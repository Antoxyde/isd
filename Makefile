DEBUG ?= 1
ifeq ($(DEBUG), 1)
	CFLAGS=-g -Wextra -Werror -Wall $(INCDIRS) --std=c99 -Wl,-rpath=. -Wno-unused-function -mbmi
else
	CFLAGS=-march=native -Ofast $(INCDIRS) --std=c99 -Wl,-rpath=. -mbmi
endif


CC=gcc
INCDIRS=-I.
LDFLAGS=
EXECUTABLES=main_stern main_prange

all: $(EXECUTABLES)

%.o : %.c %.h
	$(CC) $(CFLAGS) -c $<

main_stern.o: main_stern.c stern.h utils.h libpopcnt.h xoshiro256starstar.h
main_prange.o: main_prange.c prange.h utils.h libpopcnt.h xoshiro256starstar.h
prange.o: prange.c prange.h utils.h libpopcnt.h xoshiro256starstar.h
stern.o: stern.c utils.h libpopcnt.h xoshiro256starstar.h

main_stern: utils.o main_stern.o libm4ri-0.0.20200125.so stern.o
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

main_prange: utils.o main_prange.o prange.o libm4ri-0.0.20200125.so
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)


clean:
	rm -f $(EXECUTABLES) *.o

