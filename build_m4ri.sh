#!/bin/sh

cd m4ri
autoreconf --install
./configure
make
make check
cp .libs/libm4ri-*.so ../libm4ri.so
