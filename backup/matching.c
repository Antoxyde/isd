#include <stdio.h>

int main(void) {

    int t1[] = {1,2,3,4,4,4,4,5,5, 5, 5};
    int t2[] = {1,1,1,2,3,3,3,4,5, 5};

    int t1size = 11, t2size = 10;
    int it1 = 0, it2 = 0, sit2 = 0;

    for (it1 = 0; it1 < t1size; it1++) {

        if (it1 > 0 && t1[it1 - 1] != t1[it1]) {
            sit2 = it2;
        }

        for (it2 = sit2; it2 < t2size && t1[it1] == t2[it2]; it2++) {
            printf("%d,%d\n", t1[it1],t2[it2]);
        }
    }

    return 0;
}




