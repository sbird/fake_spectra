#Define this to disable HDF5
#HDF5=no
#Define this to disable reading of type 2 snapshots
GREAD=no
#Change this to where you installed GadgetReader
GREADDIR=${CURDIR}/../GadgetReader
#Python include path
PYINC=`pkg-config --cflags python3`

ifeq ($(CC),cc)
    GCC:=$(shell which gcc --tty-only 2>&1)
    #Can we find gcc?
    ifeq (/gcc,$(findstring /gcc,${GCC}))
       CC = gcc
       CXX = g++
    endif #Did we find gcc?
endif #CC=cc

#Are we using gcc or icc?
ifeq (icpc,$(findstring icpc,${CXX}))
  CFLAGS +=-O3 -g -w1 -openmp -fpic -march=native
  LINK +=${CXX} -O3 -openmp -march=native
else
  CFLAGS +=-O3 -g -Wall -fPIC -ffast-math
  LINK +=${CXX} -g -O3 $(PRO) -ffast-math
  LFLAGS += -lm -lgomp
endif
GCCV:=$(shell gcc --version )
ifneq (darwin,$(findstring darwin,${GCCV}))
  LFLAGS += -lgomp
  LINK += -fopenmp
  CFLAGS += -fopenmp
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
CFLAGS += -DNOGREAD
LINK += -DNOGREAD
else
IGREAD = -I${GREADDIR}
LIBS+=-lrgad -L${GREADDIR} -Wl,-rpath,${GREADDIR}
endif

# Voigt profiles vs. Gaussian profiles
CFLAGS += -DVOIGT
#If defined, looks for NHEP and NHEPP instead of NE
#CFLAGS += -DSPLIT_NE
#Use a top hat kernel instead of an SPH kernel
#CFLAGS += -DTOP_HAT_KERNEL
COM_INC = index_table.h absorption.h part_int.h
CXXFLAGS += $(CFLAGS)

.PHONY: all clean dist test python

all: _spectra_priv.so
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
	$(CXX) $(CXXFLAGS) ${IGREAD} -c $<

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

