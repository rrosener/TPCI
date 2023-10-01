# Makefile for compiling Cloudy on a UNIX platform
# This works only with GNU make!
# you need g++ or pathCC to create the dependencies (set with CXXDEPEND)
# but after that you can use any compiler you want (set with CXX)

#OPT = 
OPT = -O3 -ftrapping-math -fno-math-errno 
DEBUGOPT = -finline
CXX = g++
# -W is the old name for -Wextra, g++ 3.3 only supports the former
CXXFLAGS = -ansi ${OPT} -Wall -W #-g
LDFLAGS = ${OPT} -Wall #-g
LDLIBS =
NO_TRAP_FLAGS =
# this is intended for supplying extra command line parameters, e.g. extra -D flags
EXTRA =
SRCDIR = .
OBJEXT = o
LIBEXT = a

.PHONY: all default data debug clean distclean depend debug-test test valgrind-test

# Default rule comes first -- having this implement 'make all' is standard
all : default

# include a configuration file if one exists (or carry on without
# worrying if it doesn't)
-include Makefile.conf

main := maincl.cpp
fullsource := ${notdir ${wildcard ${SRCDIR}/*.cpp}}
source := ${filter-out Test%,${fullsource}}
testsource := ${filter Test%,${fullsource}}
libsource := ${filter-out ${main},${source}}
includes := ${notdir ${wildcard ${SRCDIR}/*.h}}
fulldepends := ${fullsource:.cpp=.d}
objects := ${source:.cpp=.${OBJEXT}}
libobjects := ${libsource:.cpp=.${OBJEXT}}
testobjects := ${testsource:.cpp=.${OBJEXT}}

ifneq (${SRCDIR},.)
  vpath %.cpp ${SRCDIR}
  vpath %.h   ${SRCDIR}
endif

# this is the compiler used to create the dependencies
CXXDEPEND = g++

# check compiler capabilities
CAPABILITIES := $(shell $(SRCDIR)/capabilities.pl $(CXX))
PRECOMPILE := $(filter precompile,$(CAPABILITIES))
VECTORIZE := $(filter vectorize,$(CAPABILITIES))

# DEP_GOALS will be an empty string for targets that don't need dependencies
NODEP_GOALS = clean distclean depend data
DEP_GOALS = ${filter-out ${NODEP_GOALS},${MAKECMDGOALS}}

# the default target does need dependencies
ifeq (${MAKECMDGOALS},)
  DEP_GOALS = need_dependencies
endif

ifeq ($(VECTORIZE),vectorize)
  OPT += -ftree-vectorize
endif

DEBUG_GOALS = $(filter debug%,$(MAKECMDGOALS))
ifneq ($(DEBUG_GOALS),)
  OPT = $(DEBUGOPT)
endif

CXXDEFINES = ${filter -D%,${CXXFLAGS}}

# Precompiling the cddefines.h header speeds up compilation dramatically with g++
# versions from 3.4.0 onwards -- unfortunately, only one header can be precompiled.
# See http://gcc.gnu.org/onlinedocs/gcc-4.2.1/gcc/Precompiled-Headers.html

presource := cddefines.h
substname := cddefines.h
preobjects :=
predepends :=
ifeq ($(PRECOMPILE),precompile)
  substname := cddefines.h.gch
  preobjects := ${presource:.h=.h.gch}
  predepends := ${presource:.h=.h.d}
endif

RANLIB := ${shell which ranlib}
ifneq ($(RANLIB),)
  RANLIB = ranlib
else
  RANLIB = ar ts
endif

PWD := ${shell pwd}

CLDCONFIG = $(PWD)/cloudyconfig.h

CDP = ${CLOUDY_DATA_PATH}
ifeq ($(CDP),)
# create reasonable default when CLOUDY_DATA_PATH not set...
  CDP = $(PWD)/${SRCDIR}/../data/
endif
ifneq (${CLOUDY_LAPACK_PATH},)
  CXXFLAGS += -DLAPACK -I${CLOUDY_LAPACK_PATH}/include
  LDLIBS += -L${CLOUDY_LAPACK_PATH}/lib -llapack
endif

# this is a collection of compiler flags that enable FPE traps
TRAP_FPE_FILTER = -fnonstd -ftrap=common -Ktrap=fp

# we dont need the -DSVN_REVISION here, disable FPE traps for unit tests
CXXUNITTESTFLAGS := $(filter-out ${TRAP_FPE_FILTER},${CXXFLAGS}) $(NO_TRAP_FLAGS)

SVNVER := ${shell which svnversion 2> /dev/null}
SED := ${shell which sed 2> /dev/null}
REVISION := exported
ifneq ($(SVNVER),)
  ifneq ($(SED),)
    REVISION := ${shell $(SVNVER) $(SRCDIR) | $(SED) 's/ /_/g'}
  endif
endif
CXXFLAGS += -DSVN_REVISION=\"$(REVISION)\"

CXXFLAGS += -DSYS_CONFIG=\"$(CLDCONFIG)\"

# add in user-defined extra flags
CXXFLAGS += $(EXTRA)
LDFLAGS += $(EXTRA)

# this is a filter to erase optimization flags from CXXFLAGS
OPTFILTER = -O -O1 -O2 -O3 -O4 -O5 -Os -xO -xO1 -xO2 -xO3 -xO4 -xO5 /O /O1 /O2 /O3 -fast -xvector%
CFGFILTER = -DSYS_CONFIG%

CXXFLAGSNOOPT := ${filter-out ${OPTFILTER},${CXXFLAGS}}
CXXFLAGSNOCFG := ${filter-out ${CFGFILTER},${CXXFLAGS}}

# Makefile.targets can override or add to the DEFAULT variable.
# the target "depend" needs to be included to assure that the build fails
# if creating one of the dependency files fails...
DEFAULT = depend cloudy.exe data

# include additional targets from file if one exists (or carry on without
# worrying if it doesn't)
-include Makefile.targets

debug : default

default : $(DEFAULT)

cloudy.exe : ${preobjects} maincl.${OBJEXT} libcloudy.${LIBEXT}
	${CXX} ${LDFLAGS} -o cloudy.exe maincl.${OBJEXT} -L. -lcloudy ${LDLIBS}

libcloudy.${LIBEXT}: ${libobjects}
	ar ru libcloudy.${LIBEXT} $^
	${RANLIB} libcloudy.${LIBEXT}

data:
	cd ${SRCDIR}/../data/cdms+jpl && $(MAKE)

clean :
	rm -f *.${OBJEXT}
	rm -f *.h.gch
	rm -rf SunWS_cache
	rm -f libcloudy.${LIBEXT}
	rm -f $(DEFAULT)
	rm -f runtests.exe
	rm -f tmp_cloudyconfig.*

distclean :
	rm -f *.${OBJEXT}
	rm -f *.h.gch
	rm -f *.d
	rm -f cloudyconfig.h
	rm -rf SunWS_cache
	rm -rf cloudy.exe.dSYM
	rm -f libcloudy.${LIBEXT}
	rm -f $(DEFAULT)
	rm -f runtests.exe
	rm -f tmp_cloudyconfig.*
	rm -f ${SRCDIR}/Makefile.conf
	rm -rf lib
	rm -rf include
	cd ${SRCDIR}/../library/UnitTest++ && $(MAKE) clean
	cd ${SRCDIR}/../data/cdms+jpl && $(MAKE) clean

depend: ${predepends} ${fulldepends}

# Ensure version update if *any* source files change
version.${OBJEXT}: ${main} ${libsource} ${includes}

cpu.${OBJEXT}: 
	${CXX} ${CXXFLAGS} -c -DCLOUDY_DATA_PATH=\"${CDP}\" $<

# an explicit -O0 is needed for some compilers like icc
parse_crashdo.${OBJEXT}:
	${CXX} ${CXXFLAGSNOOPT} ${DEBUGOPT} -c $<

# this default rule is needed on Windows
%.${OBJEXT}: %.cpp
	${CXX} ${CXXFLAGS} -c -o $@ $<

%.h.gch: %.h
	${CXX} ${CXXFLAGS} -o $@ $<

$(CLDCONFIG): ${SRCDIR}/configure.sh
	@echo "Creating cloudyconfig.h"
	@${SRCDIR}/configure.sh ${CXX} ${CXXFLAGSNOCFG} -I${SRCDIR}

debug-test : test

test: runtests.exe
	./runtests.exe

valgrind-test: runtests.exe
	valgrind --leak-check=full --track-fds=yes ./runtests.exe

runtests.exe: lib/libUnitTest++.${LIBEXT} libcloudy.${LIBEXT} ${testobjects}
	${CXX} ${CXXFLAGS} -${OBJEXT} $@ ${filter-out %.${LIBEXT}, $^} ${LDFLAGS} -L. -lcloudy -Llib -lUnitTest++ ${LDLIBS}

# the dependency on libUnitTest++.a is needed to make sure header files are in place.
# each sys_xxx dir has its own set of header files, this assures that if we do
# "make distclean" in another sys_xxx, it does not break the installation here...
${testobjects}: lib/libUnitTest++.${LIBEXT}
	${CXX} ${CXXUNITTESTFLAGS} -Iinclude/UnitTest++ -c -o $@ ${SRCDIR}/${@:.o=.cpp}

# the dependency is needed so that in a parallel build we are sure the library is there before we copy it
# we use TestUnitTest++ to make sure that the unit tests of the library passed as well...
lib/libUnitTest++.${LIBEXT}: lib/TestUnitTest++
	cd ${SRCDIR}/../library/UnitTest++; \
	mkdir -p ${PWD}/lib ${PWD}/include ${PWD}/include/UnitTest++ \
	${PWD}/include/UnitTest++/Posix ${PWD}/include/UnitTest++/Win32; \
	cp src/*.h ${PWD}/include/UnitTest++; \
	cp src/Win32/*.h ${PWD}/include/UnitTest++/Win32; \
	cp src/Posix/*.h ${PWD}/include/UnitTest++/Posix; \
	cp libUnitTest++.${LIBEXT} ${PWD}/lib/; \
	cp TestUnitTest++ ${PWD}/lib/

lib/TestUnitTest++:
	cd ${SRCDIR}/../library/UnitTest++; \
	$(MAKE) clean; $(MAKE) CC=${CXX} CCFLAGS="${CXXUNITTESTFLAGS}" LDFLAGS="${LDFLAGS}" LDLIBS="${LDLIBS}"

# only include dependencies when we really need them, this prevents unnecessary builds of the dependencies
ifneq (${DEP_GOALS},)
-include ${predepends} ${fulldepends}
endif

# $(CLDCONFIG) needs to be included here to avoid creating the dependencies twice when building from scratch
%.d: %.cpp $(CLDCONFIG)
	${SRCDIR}/make_depend.pl "${CXXDEPEND} -MM -MG -DMM ${CXXDEFINES} $<" "s/${presource}/${substname}/" "s/\.o:/.d:/" $@

%.h.d: %.h $(CLDCONFIG)
	${SRCDIR}/make_depend.pl "${CXXDEPEND} -MM -MG -DMM ${CXXDEFINES} $<" "s/\.o:/.h.gch:/" "s/\.o:/.h.d:/" $@
