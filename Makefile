DEBUG ?= 1
AVX ?= 1
ifeq ($(DEBUG), 1)
	CFLAGS=-g -Wextra -Werror -Wall $(INCDIRS) --std=c99 -Wl,-rpath=. -Wno-unused-function -mbmi -DDEBUG -fsanitize=address
else
	ifeq ($(AVX), 1)
		CFLAGS=-Ofast $(INCDIRS) --std=c99 -Wl,-rpath=. -mbmi -mavx512vl -mavx512f -mavx512dq
	else
	CFLAGS=-march=native -Ofast $(INCDIRS) --std=c99 -Wl,-rpath=. -mbmi
	endif
endif

CC=gcc
LDFLAGS=-lm4ri
EXECUTABLES=main_stern main_prange

all: $(EXECUTABLES)

%.o : %.c %.h
	$(CC) $(CFLAGS) -c $<

main_stern.o: main_stern.c stern.h utils.h libpopcnt.h xoshiro256starstar.h
main_prange.o: main_prange.c prange.h utils.h libpopcnt.h xoshiro256starstar.h
main_ballcoll.o: main_ballcoll.c ballcoll.h utils.h libpopcnt.h xoshiro256starstar.h

prange.o: prange.c prange.h utils.h libpopcnt.h xoshiro256starstar.h
stern.o: stern.c utils.h libpopcnt.h xoshiro256starstar.h
ballcoll.o: ballcoll.c utils.h libpopcnt.h xoshiro256starstar.h

main_stern: utils.o main_stern.o  stern.o #libm4ri-0.0.20200125.so
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

main_prange: utils.o main_prange.o prange.o #libm4ri-0.0.20200125.so
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

clean:
	rm -f $(EXECUTABLES) *.o

nns:
	gcc -o nns_test nns_test.c -lm4ri
