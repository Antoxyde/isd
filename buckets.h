#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

typedef struct bucket_ {
    uint64_t* tab;
    uint64_t curlen;
    uint64_t maxlen;
} bucket;

bucket** bucket_init(size_t nb_buckets, size_t bucketsize, size_t bucketinc);
void bucket_put(bucket** buckets,  uint64_t key, uint64_t data);
void bucket_free(bucket** buckets, size_t nb_buckets);
void bucket_print(bucket** buckets, size_t nb_buckets);
