#include "op339.h"
#include <stdio.h>

uint64_t mul339(const uint64_t a, const uint64_t b) {

	uint64_t a1 = a >> 17ULL;
	uint64_t a0 = a & 0x1ffffULL;

	uint64_t b1 = b >> 17ULL;
	uint64_t b0 = b & 0x1ffffULL;

	uint64_t tmp = a1 * b1 * 18 + ((a0 * b1 + b0 * a1) << 17ULL) + (b0 * a0);

	uint64_t res = 9 * (tmp >> 33ULL) + (tmp & 0x1ffffffffULL);

	if (res >= 8589934583ULL) {
		res -= 8589934583ULL;
	}

	return res;

}

uint64_t hash339(const uint32_t k, const void* buf, const size_t buflen) {

	unsigned int i;
	uint32_t* m = (uint32_t*)buf;
	uint64_t res = 0;

	for (i = 0; i < buflen / 4; i++) {
		res = mul339(k , m[i] + res);
	}

	int n = buflen & 3; // modulo 4

	if (n != 0) {

		// Cas ou on a un bloc de moins de 4 octet a la fin
		// On calcule sa taille (ici n) et on le copie dans une zone mémoire
		// ou on pourra le cast en uint32_t sans risquer d'avoir un ou plusieurs
		// octets aléatoire à la fin.

		uint32_t m_last = 0;

		uint8_t* m_cast = ((uint8_t*)m);

		for (i = buflen - n; i < buflen; i++) {
			m_last ^= (((uint32_t)m_cast[i]) << (8 * (i - buflen + n)));
		}

		res = mul339(k , m_last + res);

	}

	return res;

}
void tests_mul339(void) {

	printf("Execution des tests sur mul339 ...\n");

	//uint32_t k = 42; // tests can't fail with 42 :^)

	// Test N°1 : 2 entiers dont le résultat de leur multilpication dans Z serait supérieur a 2^64
	uint64_t a = 8589934592ULL; // 2^33
	uint64_t b = 8589934592ULL; // 2^33

	// Résultat attendu :
	// >>> ((2**33) * (2**33)) % ((2**33) - 9)
	// 81
	uint64_t res = 81ULL;

	printf("Test n°1 : %s\n", mul339(a,b) == res ? "\e[0;32mOK\e[0m" : "\e[1;31mNOPE\e[0m");

	// Test N°2 : 2 entiers dont un nul.
	a = 0ULL;
	b = 12312ULL;
	res = 0ULL;
	printf("Test n°2 : %s\n", mul339(a,b) == res ? "\e[0;32mOK\e[0m" : "\e[1;31mNOPE\e[0m");

	// Test N°3 : 2 entier, un supérieur à 2^33 - 9 et l'autre valant 1
	// Résultat attendu : la valeur du premier mod 2^33 - 9
	// >>> 2**33 + 2**25 + 4
	// 8623489028
	// >>> (2**33 + 2**25 + 4) % ((2**33) - 9)
	// 33554445

	a = 1ULL;
	b = 8623489028ULL;
	res = 33554445ULL;
	printf("Test n°3 : %s\n", mul339(a,b) == res ? "\e[0;32mOK\e[0m" : "\e[1;31mNOPE\e[0m");

}

void tests_hash339(void) {

	printf("Execution des tests sur hash339 ...\n");
	uint32_t k = 42;

	// Test N°1 : 2 blocs de 4 octets
	uint8_t t[] = {0x11, 0x22, 0x33, 0x44, 0x55, 0x66,0x77, 0x88};

	// Résultat attendu :
	// >>> (((0x44 << (3*8)) + (0x33 << (2*8)) + (0x22 << 8) + 0x11) * 42**2 + ((0x88 << (3*8)) + (0x77 << (2*8)) + (0x66 << 8) + 0x55) * 42) % ((2**33) - 9)
	// 1408077756
	uint64_t res = 1408077756ULL;

	printf("Test n°1 : %s\n", hash339(k, t, 8) == res ? "\e[0;32mOK\e[0m" : "\e[1;31mNOPE\e[0m");

	// Test N°2 : un nombre d'octet non divisible par 4
	uint8_t t2[] = {0x11, 0x22, 0x33, 0x44, 0x55, 0x66};

	// Résultat attendu :
	//>>> (((0x44 << (3*8)) + (0x33 << (2*8)) + (0x22 << 8) + 0x11) * (42**2) + (((0x66 << (1*8)) + (0x55 << (0*8))) * 42)) % ((2**33) - 9)
	// 8328286032
	res = 8328286032ULL;

	printf("Test n°2 : %s\n", hash339(k, t2, 6) == res ? "\e[0;32mOK\e[0m" : "\e[1;31mNOPE\e[0m");

}


