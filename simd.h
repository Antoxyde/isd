#ifndef ISD_POPCNT_H
#define ISD_POPCNT_H


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

#endif
