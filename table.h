#ifndef _HASHTABLE_
#define _HASHTABLE_

#include <stdint.h>
#include <stddef.h> // size_t

typedef struct bucket_ {
	rci_t** data; // We store the p indexes of the linear combination
	size_t len; // number of element currently in the bucket
	size_t count; // the bucket currently have allocated memory for `count` elements
} bucket;

typedef struct table_ {
	bucket** buckets;
	size_t tablen; // Nombre de buckets
	size_t base_bucketlen; // Taille de base d'un bucket
} table;

table* table_init(const size_t tablen, const size_t base_buckletlen);
void bucket_free(bucket* b);
void bucket_free_full(bucket* b);

void table_free(table* t);
void table_free_full(table* t);

void table_insert(table* t, const void* data, const size_t datalen, const uint64_t key);
elem* table_retrieve(const table* ht, const uint64_t key);

void table_print(const table* t);
void bucket_print(const bucket* b);

#endif
