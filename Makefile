CC=gcc
INCDIRS=-Im4ri/m4ri/ -Im4ri
LDFLAGS=
override CFLAGS += -Wextra -Werror -Wall $(INCDIRS) --std=c99 -Wl,-rpath=.
EXECUTABLES=main

all: $(EXECUTABLES)

%.o : %.c %.h
	$(CC) $(CFLAGS) -c $<

main.o: main.c prange.h utils.h libpopcnt.h xoshiro256starstar.h
prange.o: prange.c prange.h utils.h libpopcnt.h xoshiro256starstar.h
hashtable.o: hashtable.c hashtable.h utils.h op339.h
stern.o: stern.h stern.c hashtable.h utils.h
main: utils.o main.o prange.o libm4ri-0.0.20200125.so op339.o hashtable.o stern.o
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

clean:
	rm -f $(EXECUTABLES) *.o

