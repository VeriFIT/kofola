# Kofola: Modular Complementation and Inclusion Checking for Omega Automata

![build workflow](https://github.com/VeriFIT/kofola/actions/workflows/build.yml/badge.svg)

Kofola is an open source tool for an efficient complementation and inclusion checking 
of automata over infinite words (omega automata). Kofola has been built on top of [SPOT](https://spot.lrde.epita.fr/) and
inspired by [Seminator](https://github.com/mklokocka/seminator) and
[COLA](https://github.com/liyong31/COLA). 


### Requirements and dependencies

For a successful build of Kofola, `cmake` of version 3.16 (or higher) together with a C++ compiler with a support of C++-17 standard is required. Additional requirements 
include:

* [Spot](https://spot.lrde.epita.fr/)

For Debian-like systems, Spot can be installed using the package `libspot-dev`. 
For building Spot from sources, download the recent version and run:  
```
./configure
```
One can set the maximal number of colors for an automaton when configuring Spot with `--enable-max-accsets=INT`.
```
make
sudo make install
```

### Building Kofola
Please run the following steps to compile Kofola after cloning this repo:
```
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

Then you will get an executable in `build/src/kofola`. Alternatively you can 
run `make release` in the root directory.

