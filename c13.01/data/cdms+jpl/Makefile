# Makefile for converting the CDMS / JPL data files
# use 'override' here to make sure we do not inherit the settings from the main compilation
# these could be for an incompatible compiler...
override CXX = g++
override OPT = -O3 -ftrapping-math -fno-math-errno
# -W is the old name for -Wextra, g++ 3.3 only supports the former
override CXXFLAGS = -ansi ${OPT} -Wall -W -g
override EXTRA =

.PHONY: all clean squeaky-clean

data := ${notdir ${wildcard *.cat}}
processed := ${data:.cat=.dat}

all: ${processed}

convert_calpgm: convert_calpgm.cpp
	${CXX} ${CXXFLAGS} ${EXTRA} -o convert_calpgm convert_calpgm.cpp

%.dat: %.cat convert_calpgm
	./convert_calpgm $< > $@

clean:
	rm -f convert_calpgm
	rm -rf convert_calpgm.dSYM
	rm -f c??????.dat

squeaky-clean:
	rm -f convert_calpgm
	rm -f c??????.dat
	rm -f partition_function_cdms.html
	rm -f partition_function_jpl.html
	rm -f c??????.cat
	rm -f masterlist
