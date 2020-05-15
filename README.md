# ISD


# TODO

* [X] Naive Prange/Lee-Brickell implementation
* [ ] Naive Stern implementation
* [ ] SIMD-ed popcnt / linear combination
* [X] Naive popcnt without AVX512 so you can run it regardless of your CPU
* [ ] Proper makefile with AVX512 flag (atm to compile with avx512, you have to use `make "CFLAGS=-mavx512f -mavx512bw"` which is kinda ugly)

# Dependencies

* [m4ri](https://bitbucket.org/malb/m4ri/src/master/)


