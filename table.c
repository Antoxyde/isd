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
	if (t) {

		for (i = 0; i < t->nb_buckets; i++) {
			if (t->buckets[i]) {
				bucket_free_full(t->buckets[i]);
			}
		}
		free(t->buckets);
		free(t);
	}
}

void table_insert(table* t, const int* data, const size_t datalen, const uint64_t key) {

    size_t i;
	bucket* b = t->buckets[key % t->nb_buckets];
	elem* e = (elem*) malloc(sizeof(elem));

	e->data = malloc(datalen);
	CHECK_MALLOC(e->data);
	memcpy(e->data, data, datalen * sizeof(int));
	e->datalen = datalen;
    e->key = key;

	if (b == NULL) {

		b = (bucket*)malloc(sizeof(bucket));
		CHECK_MALLOC(b);

		b->len = t->base_bucketlen;
		b->count = 1;
		b->elems = (elem**)calloc(t->base_bucketlen, sizeof(elem*));
		CHECK_MALLOC(b->elems);
		b->elems[0] = e;

		t->buckets[key % t->nb_buckets] = b;
        printf("Now t->buckeys[%lu] = %p\n", key % t->nb_buckets, b);


	} else {

		if (b->count < b->len) {

			for (i = 0; i < b->len && b->elems[i]; i++);

			b->elems[i] = e;
			b->count++;

		} else {

			// No more place in the bucket, we have ro realloc
			i = b->len;

			b->len += t->base_bucketlen;

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
		for (i = 0; i < t->nb_buckets; i++) {
			if (t->buckets[i]) {
				printf("Bucket [%ld] => {\n", i);
				bucket_print(t->buckets[i]);
				printf("}\n");
			}
		}
	}
}

void bucket_print(const bucket* b) {
	size_t i;
	if (b) {
		for (i = 0; i < b->len; i++) {
			if (b->elems[i]) {
				printf("bucket @%p[%ld] => k=%lu\n", b, i, b->elems[i]->key);
			}
		}
	}
}

bucket* table_retrieve_bucket(const table* t, const uint64_t key) {
	bucket* b = t->buckets[key % t->nb_buckets];
    return b;
}

int main(void) {

    table* t = table_init(262144, 1);

    for (uint64_t i = 0; i < 20; i++) {
        table_insert(t, (int*)&i, 1, i % 10);
    }

    for (uint64_t i = 0; i < 10; i++)
        bucket_print(table_retrieve_bucket(t, i));
}



