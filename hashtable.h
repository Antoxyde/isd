#ifndef _HASHTABLE_
#define _HASHTABLE_

#include <stdint.h>
#include <stddef.h> // size_t

typedef struct elem_ {
	void* key;
	size_t keylen;
	void *data;
	size_t datalen;
	uint64_t hash_key;
} elem;

typedef struct bucket_ {
	elem** elems;
	size_t len;
	size_t count;
} bucket;

typedef struct hashtable_ {
	bucket** table;
	size_t tablen; // Nombre de buckets
	size_t base_bucketlen; // Taille de base d'un bucket
	uint32_t k;
} hashtable;


hashtable* hashtable_init(const size_t tablen, const size_t base_buckletlen);
void bucket_free(bucket* b);
void bucket_free_full(bucket* b);

void hashtable_free(hashtable* ht);
void hashtable_free_full(hashtable* ht);

void hashtable_insert(hashtable* ht, const void* data, const size_t datalen, const void* key, const size_t keylen);
elem* hashtable_retrieve(const hashtable* ht, const void* key, const size_t keylen);
int hashtable_delete(hashtable* ht, const void* key, const size_t keylen);

void hashtable_print(const hashtable* ht);
void bucket_print(const bucket* b);

#endif
