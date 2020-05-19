CC=gcc
INCDIRS=-Im4ri/m4ri/ -Im4ri
LDFLAGS=
override CFLAGS += -Wextra -Werror -Wall $(INCDIRS) --std=c99 -O2 -Wl,-rpath=.
EXECUTABLES=main

all: $(EXECUTABLES)

%.o : %.c %.h
	$(CC) $(CFLAGS) -c $<

main.o: main.c isd.h utils.h libpopcnt.h
isd.o: isd.c isd.h utils.h libpopcnt.h
main: utils.o main.o isd.o libm4ri-0.0.20200125.so
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

clean:
	rm -f $(EXECUTABLES) *.o

