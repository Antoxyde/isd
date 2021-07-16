EXECUTABLES=main_stern
CC=gcc
CFLAGS=--std=c99 -mbmi

all: 
	echo -e "Either run : \nmake local : release mode, without avx\nmake debug: debug mode, without avx\nmake avx: release mode, with avx"

avx: CFLAGS +=  -mbmi -Ofast -DFILTER -mavx512vl -mavx512f -mavx512dq -Im4ri
avx: LDFLAGS += '-Wl,-rpath,$$ORIGIN' -L. -l:libm4ri-0.0.20200125.so
avx: $(EXECUTABLES)

local: CFLAGS += -mbmi -Ofast
local: LDFLAGS += -lm4ri
local: $(EXECUTABLES)

debug: CFLAGS += -g -Wextra -Wall -Wno-unused-function -DDEBUG -fsanitize=address
debug: LDFLAGS += -lasan -lm4ri
debug: $(EXECUTABLES)

nns_test.o: nns_test.c utils.h libpopcnt.h 
main_stern.o: main_stern.c stern.h utils.h libpopcnt.h xoshiro256starstar.h
main_prange.o: main_prange.c prange.h utils.h libpopcnt.h xoshiro256starstar.h

prange.o: prange.c prange.h utils.h libpopcnt.h xoshiro256starstar.h
stern.o: stern.c utils.h libpopcnt.h xoshiro256starstar.h
semi_bc.o: semi_bc.c utils.h libpopcnt.h xoshiro256starstar.h

%.o : %.c %.h
	$(CC) $(CFLAGS) -c $<

main_stern: utils.o main_stern.o  stern.o combinations.o
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)
   
main_semibc: utils.o main_semibc.o  semi_bc.o 
	$(CC) -o $@ $^ $(FLAGS) $(LDFLAGS)

main_prange: main_prange.o prange.o  utils.o 
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

nns_test: utils.o nns_test.o  buckets.o
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)
clean:
	rm -f $(EXECUTABLES) *.o

