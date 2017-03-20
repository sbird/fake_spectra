#Optional file formats for the C reader.
#Does not affect python.
#Define this to disable HDF5
#HDF5=no
#Define this to disable reading of type 2 snapshots
GREAD=no
#Change this to where you installed GadgetReader
GREADDIR=${CURDIR}/../GadgetReader

# Voigt profiles vs. Gaussian profiles
CFLAGS += -DVOIGT

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
  LFLAGS += -lm
endif
GCCV:=$(shell gcc --version )
ifneq (darwin,$(findstring darwin,${GCCV}))
  LFLAGS += -lgomp
  LINK += -fopenmp
  CFLAGS += -fopenmp
endif
#Mac's ld doesn't support --no-add-needed, so check for it.
#We are looking for the response: ld unknown option: --no-add-needed
LDCHECK:=$(shell ld --as-needed 2>&1)
ifneq (unknown,$(findstring unknown,${LDCHECK}))
  PYLIB +=-Wl,--no-add-needed,--as-needed
endif
#Set up the python paths, if people are using non-standard python
#installs such as anaconda (very common on mac)
#If you want to use python2, change the below to python2-config.
PYPREF=$(shell python3-config --prefix)
PYLIB +=$(shell python3-config --ldflags)
#Python include path
PYINC:=-I${PYPREF}/include $(shell python3-config --includes)
#Set include path for numpy
PYINC+= -I$(shell python3 _np_setup.py)

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

COM_INC = index_table.h absorption.h part_int.h
CXXFLAGS += $(CFLAGS)

.PHONY: all clean dist test python

all: _spectra_priv.so
python: _spectra_priv.so

extract: cextract/main.o cextract/read_snapshot.o cextract/read_hdf_snapshot.o absorption.o cextract/init.o index_table.o part_int.o Faddeeva.o
	$(LINK) $(LFLAGS) $(LIBS) $^ -o $@

statistic: cextract/statistic.o cextract/calc_power.o cextract/mean_flux.o cextract/smooth.o cextract/powerspectrum.o
	$(LINK) $(LFLAGS) -lfftw3 $^ -o $@

cextract/%.o: cextract/%.c cextract/global_vars.h
cextract/main.o: cextract/main.cpp cextract/global_vars.h $(COM_INC)

cextract/read_snapshot.o: cextract/read_snapshot.cpp cextract/global_vars.h
	$(CXX) $(CXXFLAGS) ${IGREAD} -c $< -o $@

%.o: %.cpp %.h

absorption.o: absorption.cpp absorption.h singleabs.h

part_int.o: part_int.cpp part_int.h absorption.h index_table.h

calc_power.o: calc_power.c smooth.o powerspectrum.o

py_module.o: py_module.cpp $(COM_INC)
	$(CXX) $(CFLAGS) -fno-strict-aliasing -DNDEBUG $(PYINC) -c $< -o $@

_spectra_priv.so: py_module.o absorption.o index_table.o part_int.o Faddeeva.o
	$(LINK) $(LFLAGS) $(PYLIB) -shared $^ -o $@

clean:
	rm -f *.o *.pyc extract statistic _spectra_priv.so cextract/*.o

dist: Makefile
	tar -czf flux_extract.tar.gz *.c *.h *.cpp *.py $^

btest: test.cpp absorption.o index_table.o Faddeeva.o
	${LINK} ${LFLAGS} $^ -lboost_unit_test_framework -o $@

test: btest
	@./btest

