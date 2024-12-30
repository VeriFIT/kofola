# Kofola: Modular Complementation and Inclusion Checking for Omega Automata

![build workflow](https://github.com/VeriFIT/kofola/actions/workflows/build.yml/badge.svg)

Kofola is an open source tool for an efficient complementation and inclusion checking 
of automata over infinite words (omega automata). Kofola has been built on top of [SPOT](https://spot.lrde.epita.fr/) and
inspired by [Seminator](https://github.com/mklokocka/seminator) and
[COLA](https://github.com/liyong31/COLA). 


## Requirements and dependencies

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

## Building Kofola
Please run the following steps to compile Kofola after cloning this repo:
```
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

Then you will get an executable in `build/src/kofola`. Alternatively you can 
run `make release` in the root directory.

## Basic usage
Kofola assumes input omega automata in HOA format. The following command 
translates general (nondeterministic) omega-automaton stored in file `A.hoa` into a complementary 
omega automaton and prints it to the standard output:

```
./kofola A.hoa --complement
```

Note that Kofola produces automata with a general accepting condition (different output type might be specified 
by `--tba` for transition-based Büchi automata or `--tgba` for transition-based 
generalized Büchi automata). 

The following command then checks if the language specified by the omega automaton `A.hoa`
is included in the language specified by `B.hoa` and prints the result to the standard output:

```
./kofola A.hoa B.hoa --inclusion
```

Additional parameters might be passed using `--params`, e.g., `--params=merge_iwa=True`. 
In order to get a program help, run

```
./kofola --help
```

## Parameters

The complementation and the inclusion checking might be adjusted by the following key-value parameters passed as `--params='key1=val1;key2=val2'`:

| Key         | Value           | Description     |
| :---        | :---            | :---            |
| `merge_iwa` | `yes`,`no`  | Merge inherently weak components for the synchronous construction |
| `merge_det` | `yes`,`no`  | Merge deterministic components for the synchronous construction |
| `preproc_incl_A` | `low`,`high`  | Level of preprocessing applied on the first automaton (`--inclusion` only) |
| `preproc_incl_B` | `low`,`high`  | Level of preprocessing applied on the second automaton (`--inclusion` only) |
| `nac-alg` | `subs_tup`,`rank`  | Algorithm applied on nondeterministic accepting components.`subs_tup`=subset tuple construction, `rank`=rank-based complementation (experimental), default=determinization-based construction |
| `postponed` | `yes`,`no`  | Use postponed construction instead the synchronous one (`--complement` only) |
| `dir_sim` | `yes`,`no`  | Use direct simulation for macrostate pruning |

## Publications
- V. Havlena, O. Lengál, Y. Li, B. Šmahlíková and A. Turrini. [Modular Mix-and-Match Complementation of Büchi Automata](https://link.springer.com/chapter/10.1007/978-3-031-30823-9_13). In *Proc. of TACAS'23*, volume 13993 of LNCS, pages 249-270, 2023. Springer. 