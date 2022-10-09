## Copyright (C) 2019-2020  The Seminator Authors
## Copyright (C) 2022  The COLA Authors
##
## This file is a part of COLA, a tool for determinization
## of Buchi automata.
##
## COLA is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## COLA is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## The COLA Authors
##
## Removed some unnecessary files (e.g. the python binding for Seminator)

bin_PROGRAMS = kofola
lib_LTLIBRARIES = src/libkofola.la

AM_CPPFLAGS = -I$(srcdir)/src -I$(SPOTPREFIX)/include

kofola_LDADD = -L$(SPOTPREFIX)/lib src/libkofola.la -lspot -lbddx
src_libkofola_la_LIBADD = -L$(SPOTPREFIX)/lib -lspot -lbddx

src_libkofola_la_SOURCES =			\
  src/kofola.hpp			\
  src/kofola.cpp			\
  src/abstract_complement_alg.cpp     \
  src/complement_alg_mh.cpp     \
  src/complement_alg_ncsb.cpp     \
  src/complement_alg_safra.cpp     \
  src/complement_alg_rank.cpp     \
  src/complement_tnba.cpp				\
  src/composer.cpp				\
  src/congr.cpp				\
  src/decomposer.cpp			\
  src/determinize.cpp			\
  src/determinize_teba.cpp			\
  src/determinize_tldba.cpp			\
  src/determinize_tnba.cpp			\
  src/determinize_twba.cpp			\
  src/optimizer.cpp				\
  src/simulation.cpp			\
  src/rankings.cpp        

kofola_SOURCES = src/main.cpp
