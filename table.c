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

	bucket* b = h->buckets[key];
	elem* e = (elem*) malloc(sizeof(elem));

	e->data = malloc(datalen);
	CHECK_MALLOC(e->data);
	memcpy(e->data, data, datalen);
	e->datalen = datalen;

	if (b == NULL) {

		b = (bucket*)malloc(sizeof(bucket));
		CHECK_MALLOC(b);

		b->len = t->base_bucketlen;
		b->count = 1;
		b->elems = (elem**)calloc(t->base_bucketlen, sizeof(elem*));
		CHECK_MALLOC(b->elems);

		ht->table[key % t->tablen] = b;

		b->elems[0] = e;

	} else {

		if (b->count < b->len) {

			for (i = 0; i < b->len && b->elems[i]; i++);

			b->elems[i] = e;
			b->count++;

		} else {

			// No more place in the bucket, we have ro realloc
			i = b->len;

			b->len += ht->base_bucketlen;

			b->elems = (elem**)realloc(b->elems, b->len * sizeof(elem*));
			CHECK_MALLOC(b->elems);

			for (; i < b->len; i++) {
				b->elems[i] = NULL;
			}

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

