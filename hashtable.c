#include "hashtable.h"
#include "op339.h"
#include "utils.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


hashtable* hashtable_init(const size_t tablen, const size_t base_bucketlen) {

	hashtable* ht = (hashtable*)malloc(sizeof(hashtable));
	ht->base_bucketlen = base_bucketlen;
	ht->tablen = tablen;

	do {
		ht->k = rand();
	} while (ht->k == 0 || ht->k == 1);

	ht->table = (bucket**)calloc(tablen , sizeof(bucket*)); // Calloc donc tout les pointeurs sont init a NULL

	return ht;
}

void bucket_free(bucket* b) {

	if (b) {
		free(b);
	}

}

void bucket_free_full(bucket* b) {
	size_t i;

	if (b) {
		for (i = 0; i < b->count; i++) {
			if (b->elems[i]) {
				free(b->elems[i]->data);
				free(b->elems[i]->key);
				free(b->elems[i]);
			}
		}

		free(b->elems);
		free(b);

	}
}

void hashtable_free(hashtable* ht) {

	if (ht) {
		free(ht);
	}
}

void hashtable_free_full(hashtable* ht) {
	size_t i;

	if (ht) {

		for (i = 0; i < ht->tablen; i++) {
			if (ht->table[i]) {
				bucket_free_full(ht->table[i]);
			}
		}
		free(ht->table);
		free(ht);

	}
}

void hashtable_insert(hashtable* ht, const void* data, const size_t datalen, const void* key, const size_t keylen) {

	// Ajoute l'élement ayant pour data `data` dans le bucket correspondant a la clé `key`.
	// Si l'élement est déjà présent, ne fait rien.

	// 2 cas :
	// 	- Aucun élement ayant le même hash n'a déja été inséré, on crée le bucket
	// 	- Au moins un élément est dans le bucket, on vérifie si aucun élément n'as la même clé avant d'insérer
	//	Si l'élement n'est pas déjà présent, 2 cas :
	// 		- Le nombre d'élement est entre un et base_bucketlen - 1 éléments : on a juste a inséré l'élement
	// 		- Il y a déjà base_bucketlen élement dans le bucket, on doit alloué plus de place.

    size_t i;
	uint64_t hashkey = hash339(ht->k, key, keylen);

	// Récupération du bucket correspondant
	bucket* b = ht->table[hashkey %  ht->tablen];

	// Vérifie que l'élement n'est pas déja dans le bucket
	if (b != NULL) {
		for (i = 0; i < b->len; i++) {
			if (b->elems[i] && b->elems[i]->hash_key == hashkey) {
				return;
			}
		}
	}

	// Il n'y est pas, on le crée
	elem* e = (elem*) malloc(sizeof(elem));
	e->hash_key = hashkey;

	e->data = malloc(datalen);
	CHECK_MALLOC(e->data);
	memcpy(e->data, data, datalen);
	e->datalen = datalen;

	e->key = malloc(keylen);
	CHECK_MALLOC(e->key);
	memcpy(e->key, key, keylen);
	e->keylen = keylen;

	if (b == NULL) {

		// 1er cas, le bucket n'existe pas encore, on le crée
		// on lui ajoute l'élement à ajouter
		// puis on l'ajoute à notre hashtable

		b = (bucket*)malloc(sizeof(bucket));
		CHECK_MALLOC(b);

		b->len = ht->base_bucketlen;
		b->count = 1;
		b->elems = (elem**)calloc(ht->base_bucketlen, sizeof(elem*));
		CHECK_MALLOC(b->elems);

		ht->table[hashkey % ht->tablen] = b;

		b->elems[0] = e;

	} else {

		// 2nd cas, le bucket existe déjà

		if (b->count < b->len) { // On a encore de la place dans le bucket, on ajoute simplement l'élement

			// Récupération de l'indice du premier élement NULL dans notre bucket
			for (i = 0; i < b->len && b->elems[i]; i++);

			// Place l'élement a cet indice là
			b->elems[i] = e;
			b->count++;

		} else {

			// On a plus de place dans ce bucket , il faut en reallouer
			i = b->len;

			b->len += ht->base_bucketlen;

			b->elems = (elem**)realloc(b->elems, b->len * sizeof(elem*));
			CHECK_MALLOC(b->elems);

			// Comme realloc n'initalise pas la mémoire à 0, il faut le faire manuellement
			for (; i < b->len; i++) {
				b->elems[i] = NULL;
			}

			// Place l'élement a la suite, au début du nouveau tableau reallouer
			b->elems[b->count] = e;
			b->count++;

		}
	}

}

void hashtable_print(const hashtable* ht) {

	// Affiche tout les buckets de la hashtable `ht`.

	size_t i;

	if (ht) {
		printf("[PrintHashtable] Buckets non nuls : \n");
		for (i = 0; i < ht->tablen; i++) {
			if (ht->table[i]) {
				printf("Bucket [%ld] => {\n", i);
				bucket_print(ht->table[i]);
				printf("}\n");
			}
		}
	}
}

void bucket_print(const bucket* b) {

	// Affiche chaque élement non nul du bucket `b`.

	size_t i;

	if (b) {

		for (i = 0; i < b->len; i++) {
			if (b->elems[i]) {
				printf("\t[%ld] => %lu\n", i, b->elems[i]->hash_key);
			}
		}
	}
}

elem* hashtable_retrieve(const hashtable* ht, const void* key, const size_t keylen) {

	// Si un élement de `ht` possède la clé `key`, le retourne.
	// Sinon, renvoie NULL.

	// Récupération de la hashkey liée à la clé `key`
	uint64_t hashkey = hash339(ht->k, key, keylen);

	// Récupération du bucket associé
	bucket* b = ht->table[hashkey % ht->tablen];

	if (b) { // Si le bucket est non-vide
		size_t i;

		for (i = 0; i < b->len; i++) { // On regarde chaque élement
			if (b->elems[i] && b->elems[i]->hash_key == hashkey) {
				return b->elems[i];
			}
		}
	}

	return NULL;

}

int hashtable_delete(hashtable* ht, const void* key, const size_t keylen) {

	// Si un élement de `ht` à la clé `key`, le supprime de la hashtable `ht` (free + mise a NULL du pointeur) et renvoie 0
	// Sinon, ne fait rien et renvoie 1.

	uint64_t hashkey = hash339(ht->k, key, keylen);
	bucket* b = ht->table[hashkey % ht->tablen];

	if (b) {
		size_t i;

		for (i = 0; i < b->len; i++) {
			if (b->elems[i] && b->elems[i]->hash_key == hashkey) {
				free(b->elems[i]->data);
				free(b->elems[i]->key);
				free(b->elems[i]);
				b->elems[i] = NULL;
				return 0;
			}
		}
	}

	return 1;
}
