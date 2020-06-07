#ifndef _HASHTABLE_
#define _HASHTABLE_

#include <stdint.h>
#include <stddef.h> // size_t

#define CHECK_MALLOC(X) do { \
        if ((X) == NULL) {\
	        fprintf(stderr, "Malloc failed at line %d in func %s\n", __LINE__, __func__);\
	        exit(1);\
        } \
    } while (0)

typedef struct elem_ {
    uint64_t key;
    int* data; // rci_t is a typedef for int
    size_t datalen;
} elem;

typedef struct bucket_ {
	elem** elems;
	size_t len; // number of element currently in the bucket
	size_t count; // the bucket currently have allocated memory for `count` elements
} bucket;

typedef struct table_ {
	bucket** buckets;
	size_t nb_buckets; // Nombre de buckets
	size_t base_bucketlen; // Taille de base d'un bucket
} table;

table* table_init(const size_t tablen, const size_t base_buckletlen);
void bucket_free(bucket* b);
void bucket_free_full(bucket* b);

// Free all the table elements and reset all bucket to NULL
void table_reset(table* t);

void table_free(table* t);
void table_free_full(table* t);

void table_insert(table* t, const int* data, const size_t datalen, const uint64_t key);
bucket* table_retrieve_bucket(const table* t, const uint64_t key);

void table_print(const table* t);
void bucket_print(const bucket* b);

#endif
