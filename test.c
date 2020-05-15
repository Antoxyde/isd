#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "simd.h"
#include "xoshiro256starstar.h"


int main(void) {


    uint64_t* data = malloc(4*sizeof(uint64_t*));
    uint64_t sum = 0;

    for (uint64_t i = 0; i < 100000000; i++) {

        data[0] = xoshiro256starstar_random();
        data[1] = xoshiro256starstar_random();
        data[2] = xoshiro256starstar_random();
        data[3] = xoshiro256starstar_random();

        sum += popcnt(data, 4*sizeof(uint64_t));
    }

    printf("sum : %ld\n", sum);


    return 0;
}

/*

Avec -O1 : popcnt est inliné ? Fait des appels a __popcountdi2 ? O_o
real    0m2.644s
user    0m2.643s
sys     0m0.001s

Avec -O2 : popcnt inliné pareil + gcc unroll la boucle a 4 appel par tour
real    0m1.668s
user    0m1.666s
sys     0m0.002s


Avec -O3 : popcnt inliné comme -O2 , l'appel a xoshiro256starstar_random est aussi inliné
real    0m1.024s
user    0m1.023s
sys     0m0.001s

>>> ((100000000 * 32) / (1024 * 1024 * 1024)) / 1.024
2.9103830456733704

=> version naive, 2.91 Go/s

for (uint64_t i = 0; i < 100000000; i++) {

    data[0] = xoshiro256starstar_random();
    data[1] = xoshiro256starstar_random();
    data[2] = xoshiro256starstar_random();
    data[3] = xoshiro256starstar_random();

    sum += popcnt(data, 4*sizeof(uint64_t));
}

printf("sum : %ld\n", sum);

*/


