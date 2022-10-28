# Kofola: modular complementation and inclusion checking for omega automata

Kofola has been built on top of [SPOT](https://spot.lrde.epita.fr/) and
inspired by [Seminator](https://github.com/mklokocka/seminator) and
[COLA](https://github.com/liyong31/COLA).


### Requirements
* [Spot](https://spot.lrde.epita.fr/)

```
./configure --enable-max-accsets=128
```
One can set the maximal number of colors for an automaton when configuring Spot with --enable-max-accsets=INT
```
make && make install
```

### Compilation
Please run the following steps to compile Kofola after cloning this repo:
```
mkdir build && cd build
cmake ..
make
```

Then you will get an executable in `build/src/kofola`.
