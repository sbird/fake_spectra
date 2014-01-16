#Define this to disable HDF5
#HDF5=no
#Define this to disable reading of type 2 snapshots
GREAD=no
#Change this to where you installed GadgetReader
GREAD=${CURDIR}/../GadgetReader
#Python include path
PYINC=-I/usr/include/python2.6 -I/usr/include/python2.6

GCCV:=$(shell gcc --version | head -1)
ifeq (4.8,$(findstring 4.8,${GCCV}))
	CC = gcc
	CXX = g++
endif
ifeq ($(CC),cc)
  ICC:=$(shell which icc --tty-only 2>&1)
  #Can we find icc?
  ifeq (/icc,$(findstring /icc,${ICC}))
     CC = icc -vec_report0
     CXX = icpc
  else
     GCC:=$(shell which gcc --tty-only 2>&1)
     #Can we find gcc?
     ifeq (/gcc,$(findstring /gcc,${GCC}))
        CC = gcc
        CXX = g++
     endif
  endif
endif

#Are we using gcc or icc?
ifeq (icc,$(findstring icc,${CC}))
  CFLAGS +=-O3 -g -w1 -openmp -fpic -march=native
  LINK +=${CXX} -O3 -openmp -march=native
else
  CFLAGS +=-O3 -g -Wall -fopenmp -fPIC -ffast-math
  LINK +=${CXX} -g -O3 -fopenmp $(PRO) -ffast-math
  LFLAGS += -lm -lgomp
endif
PG =
LIBS=

ifeq ($(HDF5),no)
#Support for loading HDF5 files
CFLAGS += -DNOHDF5
else
LIBS+=-lhdf5 -lhdf5_hl
endif
ifeq ($(GREAD),no)
#Support for loading HDF5 files
CXXFLAGS += -DNOGREAD
else
LIBS+=-lrgad -L${GREAD} -Wl,-rpath,${GREAD}
endif

# Voigt profiles vs. Gaussian profiles
CFLAGS += -DVOIGT
#If defined, looks for NHEP and NHEPP instead of NE
#CFLAGS += -DSPLIT_NE
COM_INC = index_table.h absorption.h part_int.h
CXXFLAGS += $(CFLAGS)

.PHONY: all clean dist test python

all: extract _spectra_priv.so statistic rescale

python: _spectra_priv.so

extract: main.o read_snapshot.o read_hdf_snapshot.o absorption.o init.o index_table.o part_int.o Faddeeva.o
	$(LINK) $(LFLAGS) $(LIBS) $^ -o $@

rescale: rescale.o powerspectrum.o mean_flux.o calc_power.o smooth.o
	$(LINK) $(LFLAGS) -lfftw3 $^ -o $@

statistic: statistic.o calc_power.o mean_flux.o smooth.o powerspectrum.o
	$(LINK) $(LFLAGS) -lfftw3 $^ -o $@

%.o: %.c global_vars.h
main.o: main.cpp global_vars.h $(COM_INC)

read_snapshot.o: read_snapshot.cpp global_vars.h
	$(CXX) $(CXXFLAGS) -I${GREAD} -c $<

%.o: %.cpp %.h

absorption.o: absorption.cpp absorption.h singleabs.h

part_int.o: part_int.cpp part_int.h absorption.h index_table.h

calc_power.o: calc_power.c smooth.o powerspectrum.o

py_module.o: py_module.cpp $(COM_INC)
	$(CXX) $(CFLAGS) -fno-strict-aliasing -DNDEBUG $(PYINC) -c $< -o $@

_spectra_priv.so: py_module.o absorption.o index_table.o part_int.o Faddeeva.o
	$(LINK) $(LFLAGS) -shared $^ -o $@

clean:
	rm -f *.o  extract rescale statistic _spectra_priv.so

dist: Makefile
	tar -czf flux_extract.tar.gz *.c *.h *.cpp *.py $^

btest: test.cpp absorption.o index_table.o Faddeeva.o
	${LINK} ${LFLAGS} -lboost_unit_test_framework $^ -o $@

test: btest
	@./btest

