/*  Written in 2018 by David Blackman and Sebastiano Vigna (vigna@acm.org)

	To the extent possible under law, the author has dedicated all copyright
	and related and neighboring rights to this software to the public domain
	worldwide. This software is distributed without any warranty.

	See <http://creativecommons.org/publicdomain/zero/1.0/>. */

/* basic wrappers for convenience, PK, 2018 */

#ifndef __XOSHIRO256starstar___
#define __XOSHIRO256starstar___

#include <stdint.h>

/* This is xoshiro256** 1.0, our all-purpose, rock-solid generator. It has
   excellent (sub-ns) speed, a state (256 bits) that is large enough for
   any parallel application, and it passes all tests we are aware of.

   For generating just floating-point numbers, xoshiro256+ is even faster.

   The state must be seeded so that it is not everywhere zero. If you have
   a 64-bit seed, we suggest to seed a splitmix64 generator and use its
   output to fill s. */

static inline uint64_t __my_little_xoshiro256starstar__rotl(const uint64_t x, int k) {
	return (x << k) | (x >> (64 - k));
}


static int __my_little_init_was_done = 0;
static uint64_t __my_little_xoshiro256starstar__s[4];

/* Inits */

void __my_little_xoshiro256starstar_initialization(uint64_t iv[4])
{
	__my_little_xoshiro256starstar__s[0] = iv[0];
	__my_little_xoshiro256starstar__s[1] = iv[1];
	__my_little_xoshiro256starstar__s[2] = iv[2];
	__my_little_xoshiro256starstar__s[3] = iv[3];

	return;
}



/* This function initializes one state with a key from /dev/urandom */
void __my_little_xoshiro256starstar_unseeded_init()
{
	FILE *urd = fopen("/dev/urandom", "r");
	uint64_t iv[4] = {1,1,1,1};

	if (urd == NULL)
	{
		fprintf(stderr, "failed to initialize the little xoshiro256** prng [No file called /dev/urandom]\n");
		__my_little_xoshiro256starstar_initialization(iv);
		return;
	}

	if(1 != fread((uint8_t *)iv, 32, 1, urd))
	{
		fprintf(stderr, "failed to initialize the little8 xoshiro256** prng [Not enough p$ bytes]\n");
	}
	fclose(urd);
	__my_little_xoshiro256starstar_initialization(iv);
	return;
}


uint64_t __my_little_xoshiro256starstar__next__unsafe(void) {
	const uint64_t result_starstar = __my_little_xoshiro256starstar__rotl(__my_little_xoshiro256starstar__s[1] * 5, 7) * 9;

	const uint64_t t = __my_little_xoshiro256starstar__s[1] << 17;

	__my_little_xoshiro256starstar__s[2] ^= __my_little_xoshiro256starstar__s[0];
	__my_little_xoshiro256starstar__s[3] ^= __my_little_xoshiro256starstar__s[1];
	__my_little_xoshiro256starstar__s[1] ^= __my_little_xoshiro256starstar__s[2];
	__my_little_xoshiro256starstar__s[0] ^= __my_little_xoshiro256starstar__s[3];

	__my_little_xoshiro256starstar__s[2] ^= t;

	__my_little_xoshiro256starstar__s[3] = __my_little_xoshiro256starstar__rotl(__my_little_xoshiro256starstar__s[3], 45);

	return result_starstar;
}

uint64_t __my_little_xoshiro256starstar__next(void) {
	if (!__my_little_init_was_done)
	{
		__my_little_xoshiro256starstar_unseeded_init();
		__my_little_init_was_done = 1;
	}

	return __my_little_xoshiro256starstar__next__unsafe();
}


/* This is the jump function for the generator. It is equivalent
   to 2^128 calls to __my_little_xoshiro256starstar__next(); it can be used to generate 2^128
   non-overlapping subsequences for parallel computations. */

void __my_little_xoshiro256starstar__jump(void) {
	static const uint64_t JUMP[] = { 0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c };

	uint64_t s0 = 0;
	uint64_t s1 = 0;
	uint64_t s2 = 0;
	uint64_t s3 = 0;
	for(int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
		for(int b = 0; b < 64; b++) {
			if (JUMP[i] & UINT64_C(1) << b) {
				s0 ^= __my_little_xoshiro256starstar__s[0];
				s1 ^= __my_little_xoshiro256starstar__s[1];
				s2 ^= __my_little_xoshiro256starstar__s[2];
				s3 ^= __my_little_xoshiro256starstar__s[3];
			}
			__my_little_xoshiro256starstar__next();
		}

	__my_little_xoshiro256starstar__s[0] = s0;
	__my_little_xoshiro256starstar__s[1] = s1;
	__my_little_xoshiro256starstar__s[2] = s2;
	__my_little_xoshiro256starstar__s[3] = s3;
}



/* This is the long-jump function for the generator. It is equivalent to
   2^192 calls to __my_little_xoshiro256starstar__next(); it can be used to generate 2^64 starting points,
   from each of which jump() will generate 2^64 non-overlapping
   subsequences for parallel distributed computations. */

void __my_little_xoshiro256starstar__long_jump(void) {
	static const uint64_t LONG_JUMP[] = { 0x76e15d3efefdcbbf, 0xc5004e441c522fb3, 0x77710069854ee241, 0x39109bb02acbe635 };

	uint64_t s0 = 0;
	uint64_t s1 = 0;
	uint64_t s2 = 0;
	uint64_t s3 = 0;
	for(int i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
		for(int b = 0; b < 64; b++) {
			if (LONG_JUMP[i] & UINT64_C(1) << b) {
				s0 ^= __my_little_xoshiro256starstar__s[0];
				s1 ^= __my_little_xoshiro256starstar__s[1];
				s2 ^= __my_little_xoshiro256starstar__s[2];
				s3 ^= __my_little_xoshiro256starstar__s[3];
			}
			__my_little_xoshiro256starstar__next();
		}

	__my_little_xoshiro256starstar__s[0] = s0;
	__my_little_xoshiro256starstar__s[1] = s1;
	__my_little_xoshiro256starstar__s[2] = s2;
	__my_little_xoshiro256starstar__s[3] = s3;
}

/*
 * Aliases
 */

uint64_t xoshiro256starstar_random(void)
{
	return __my_little_xoshiro256starstar__next();
}

uint64_t xoshiro256starstar_random_unsafe(void)
{
	return __my_little_xoshiro256starstar__next__unsafe();
}

void xoshiro256starstar_random_set(uint64_t seed[4])
{
	__my_little_xoshiro256starstar_initialization(seed);
}
#endif
