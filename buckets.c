#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include "utils.h"
#include "buckets.h"

static size_t default_bucketsize = 0;
static size_t default_bucketinc = 0;

bucket** bucket_init(size_t nb_buckets, size_t bucketsize, size_t bucketinc) {

    default_bucketsize = bucketsize;
    default_bucketinc = bucketinc;

    bucket** buckets = calloc(nb_buckets,  sizeof(bucket*));
    if (!buckets) {
        fprintf(stderr, "Can't calloc %ld buckets, exiting\n", nb_buckets);
        exit(1);
    }

    return buckets;
}

void bucket_put(bucket** buckets,  uint64_t key, uint64_t data) {
    assert(key < (1ULL << 25));
    bucket* b = buckets[key];

    if (!b) {
        b = malloc(sizeof(bucket));
        b->maxlen = default_bucketsize;
        b->curlen = 0;
        b->tab = malloc(b->maxlen * sizeof(uint64_t));
        CHECK_MALLOC(b->tab);
        buckets[key] = b;
    }
    
    if (b->curlen == b->maxlen) {
        b->maxlen += default_bucketinc;
        b->tab = realloc(b->tab, b->maxlen * sizeof(uint64_t));
        CHECK_MALLOC(b->tab);
    }
    
    b->tab[b->curlen] = data;
    b->curlen++;
}

void bucket_free(bucket** buckets, size_t nb_buckets) {
    for (unsigned int i = 0; i < nb_buckets; i++) {
        if (buckets[i]) { 
            free(buckets[i]->tab);
            free(buckets[i]);
        }
    }

    free(buckets);
}

void bucket_print(bucket** buckets, size_t nb_buckets) {
    for (size_t i = 0; i < nb_buckets; i++) {
        if (buckets[i]) {
            printf("bucket %lX:\n", i);
            for (size_t j = 0; j < buckets[i]->curlen; j++) {
                printf("\t%08lX\n", buckets[i]->tab[j]);
            }
        } else {
            printf("buckets %lX is empty.\n", i);
        }
    }
}
