#include <stdio.h>
#include <stdlib.h>
#include <m4ri.h>

int main(void) {
    rci_t n = 10, i = 0;
    mzd_t* m = mzd_init(n, n);
    mzd_randomize(m);

    printf("Avant perm : \n");
    mzd_print(m);

    rci_t* perm = (rci_t*)malloc(sizeof(rci_t) * n);
    for (; i < n; i++) perm[i] = i;

    for (i = 0; i < n; i++) {
        int a = rand() % n;
        int b = rand() % n;

        mzd_row_swap(m, a, b);
        rci_t tmp = perm[a];
        perm[a] = perm[b];
        perm[b] = tmp;
    }

    printf("Milieu perm : \n");
    mzd_print(m);

    mzd_t* mm = mzd_copy(NULL, m);

    for (i = 0; i < n; i++)
        mzd_copy_row(mm, perm[i], m, i);

    printf("Après perm ; \n");
    mzd_print(mm);



    return 0;
}

