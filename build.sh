#!/bin/bash

# removing previously built binary
rm various_works

# cleaning old build
rm -rf build/

# creating build folder
mkdir build

cd build/

# preparation of Makefile
cmake ..

# compile & build
make

# copying the binary to the root project folder
cp various_works ../

# starting binary
./various_works
