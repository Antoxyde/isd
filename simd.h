#ifndef ISD_POPCNT_H
#define ISD_POPCNT_H

#include "libpopcnt.h"

static inline uint16_t my_popcnt(void* data, int size) {

    int i = 0;
    uint16_t result = 0;
    uint64_t* values_64 = (uint64_t*) data;
    int len = size / 8;
    int remaining = size % 8;

    for (i = 0; i < len; i++)
        result += __builtin_popcountll(values_64[i]);

    uint8_t* values_8 = (uint8_t*) (data + size * 8);

    for (i = 0; i < remaining; i++)
        result += __builtin_popcount(values_8[i]);

    return result;
}


/*
static inline uint64_t popcnt64(uint64_t x) {
  __asm__ ("popcnt %1, %0" : "=r" (x) : "0" (x));
  return x;
}
*/

static inline uint64_t popcnt640_a(uint64_t* data) {
    uint64_t sum =  popcnt64(data[0]);
    sum += popcnt64(data[1]);
    sum += popcnt64(data[2]);
    sum += popcnt64(data[3]);
    sum += popcnt64(data[4]);
    sum += popcnt64(data[5]);
    sum += popcnt64(data[6]);
    sum += popcnt64(data[7]);
    sum += popcnt64(data[8]);
    sum += popcnt64(data[9]);
    return sum;

}

#endif
