#include "utils.h"

int main(void) {

    mzd_t* H = mzd_init(640, 1280);
    H = load_challenge("challenges/LW_1280_0", H);
    mzd_print(H);

    return 0;
}
