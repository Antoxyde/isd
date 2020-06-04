#include "table.h"
#include "utils.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


table* table_init(const size_t nb_buckets, const size_t base_bucketlen) {

	table* t = (table*)malloc(sizeof(table));
	t->base_bucketlen = base_bucketlen;
	t->nb_buckets = nb_buckets;
	t->buckets = (bucket**)calloc(nb_buckets , sizeof(bucket*));

	return t;
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

void table_free(table* t) {
	if (t) {
		free(t);
	}
}

void table_free_full(table* t) {
	size_t i;
	if (ht) {

		for (i = 0; i < t->nb_buckets; i++) {
			if (ht->table[i]) {
				bucket_free_full(ht->table[i]);
			}
		}
		free(t->buckets);
		free(t);
	}
}

void table_insert(table* ht, const void* data, const size_t datalen, const uint64_t key) {

	// Récupération du bucket correspondant
	bucket* b = h->buckets[key];

	// Vérifie que l'élement n'est pas déja dans le bucket
	if (b != NULL) {
		for (i = 0; i < b->len; i++) {
			if (b->elems[i]) {
				return;
			}
		}
	}

	// Il n'y est pas, on le crée
	elem* e = (elem*) malloc(sizeof(elem));

	e->data = malloc(datalen);
	CHECK_MALLOC(e->data);
	memcpy(e->data, data, datalen);
	e->datalen = datalen;

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

void table_print(const table* t) {
	uint64_t i;
	if (t) {
		printf("[PrintTable] Buckets : \n");
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

elem* table_retrieve(const table* t, const uint64_t key);

	bucket* b = t->table[key % t->tablen];

	if (b) { // Si le bucket est non-vide
		size_t i;

		for (i = 0; i < b->len; i++) { // On regarde chaque élement
			if (b->elems[i] && b->elems[i]->key == key) {
				return b->elems[i];
			}
		}
	}

	return NULL;

}

int table_delete(table* t, const uint64_t key) {

	// Si un élement de `ht` à la clé `key`, le supprime de la hashtable `ht` (free + mise a NULL du pointeur) et renvoie 0
	// Sinon, ne fait rien et renvoie 1.

	bucket* b = ht->table[key];

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
