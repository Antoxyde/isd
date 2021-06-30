# ISD algorithms

Two differents algorithms are implemented :

* Prange [[Prange](#Prange)], using Canteaut en Chabaud improvements [[CC](#CC)]
* Stern [[Stern](#Stern)], also using Canteaut and Chabaud improvements [[CC](#CC)] , and some ideas of Bernstein, Lange and Peters [[BLP](#BLP)]

Both of them are implemented with AVX-512 support.

# Usage

## Compilation

3 targets : debug, run, avx.

* debug : no avx512, debug stuff on, optimizations off
* run: no avx512, optimisations on
* avx: avx512 on, optimizations on

Use `make <target>` to compile the target you want.

## 
# Credits

* [m4ri](https://bitbucket.org/malb/m4ri/src/master/)
* [libpopcnt](https://github.com/kimwalisch/libpopcnt)

# Bibliography

<a name="Prange">[Prange]</a>  Prange, E.: The use of information sets in decoding cyclic codes. IRE Transactions IT-8(1962) S5–S9

<a name="CC">[CC]</a> Canteaut, A., Chabaud, F.:  A new algorithm for finding minimum-weight words in a linear code: Application to McEliece’s cryptosystem and to narrow-sense BCH codes of length 511.  IEEE Transactions on Information Theory44(1) (January 1998) 367–378

<a name="Stern">[Stern]</a>  Jacques  Stern.   A  method  for  finding  codewords  of  small  weight. Coding theory and applications,  volume 388  of Lecture Notes in Computer Science, 1989.

<a name="BLP">[BLP]</a> Daniel J. Bernstein, Tanja Lange, Christiane Peters: [Attacking and defending the McEliece cryptosystem](https://cr.yp.to/codes/mceliece-20080807.pdf), 2008.

