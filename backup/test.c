#include <stdio.h>
#include <stdint.h>

#define STERN_GET(tab, index) ((tab[index >> 6]) >> (index & 0x3f)) & 1
#define STERN_SET_ONE(tab, index) (tab[index >> 6]) |= (1ULL << (index & 0x3f))

int main(void) {

    uint64_t tab[2] = {0,0};

    for (int i = 0; i < 128; i += 2) {
        STERN_SET_ONE(tab, i);
    }

    printf("Now tab is : %lu, %lu\n", tab[0], tab[1]);
    printf("In binary : \n0b");

    for (int i = 0; i < 128; i++) {
        if (i == 64) printf("\n0b");

        printf("%s" , STERN_GET(tab, i) ? "1" : "0");
    }
    printf("\n");


    return 0;
}


