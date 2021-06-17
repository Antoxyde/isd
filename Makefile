EXECUTABLES=nns_test main_prange main_stern
CC=gcc
CFLAGS=--std=c99 -DPROGRESS -mbmi

all: 
	echo -e "Either run : \nmake local : release mode, without avx\nmake debug: debug mode, without avx\nmake avx: release mode, with avx"

avx: CFLAGS += -Wl,-rpath=. -mbmi -DDEBUG -Ofast -DPROGRESS -mavx512vl -mavx512f -mavx512dq
avx: LDFLAGS += -llibm4ri-0.0.20200125.so
avx: $(EXECUTABLES)

run: CFLAGS += -Wl,-rpath=. -mbmi -Ofast -DPROGRESS 
run: LDFLAGS += -lm4ri
run: $(EXECUTABLES)

debug: CFLAGS += -g -fsanitize=address -Wextra -Wall -Wno-unused-function -DDEBUG
debug: LDFLAGS += -lm4ri  -lasan
debug: $(EXECUTABLES)

%.o : %.c %.h
	$(CC) $(CFLAGS) -c $<

nns_test.o: nns_test.c utils.h libpopcnt.h 
main_stern.o: main_stern.c stern.h utils.h libpopcnt.h xoshiro256starstar.h
main_prange.o: main_prange.c prange.h utils.h libpopcnt.h xoshiro256starstar.h

prange.o: prange.c prange.h utils.h libpopcnt.h xoshiro256starstar.h
stern.o: stern.c utils.h libpopcnt.h xoshiro256starstar.h

main_stern: utils.o main_stern.o  stern.o 
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

main_prange: main_prange.o prange.o  utils.o
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

nns_test: utils.o nns_test.o  buckets.o
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)
clean:
	rm -f $(EXECUTABLES) *.o

