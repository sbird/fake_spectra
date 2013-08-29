#Change this to where you installed GadgetReader
GREAD=${CURDIR}/../GadgetReader
#Python include path
PYINC=-I/usr/include/python2.6 -I/usr/include/python2.6

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
  CFLAGS +=-O3 -g -Wall -fopenmp -fPIC
  LINK +=${CXX} -g -O3 -fopenmp $(PRO)
  LFLAGS += -lm -lgomp
endif
OPTS =
PG =
# Voigt profiles vs. Gaussian profiles
OPTS += -DVOIGT
#Support for loading HDF5 files
OPTS += -DHDF5
#This is misnamed: in reality it looks for NE instead of NHEP and NHEPP
OPTS += -DGADGET3 #-DJAMIE
CFLAGS += $(OPTS)
CXXFLAGS += $(CFLAGS)
EXT_INC = global_vars.h index_table.h absorption.h part_int.h
LIBS=-lrgad -L${GREAD} -Wl,-rpath,${GREAD} -lhdf5 -lhdf5_hl

.PHONY: all clean dist test

all: extract statistic rescale	_spectra_priv.so

extract: main.o read_snapshot.o read_hdf_snapshot.o absorption.o init.o index_table.o part_int.o
	$(LINK) $(LFLAGS) $(LIBS) $^ -o $@

rescale: rescale.o powerspectrum.o mean_flux.o calc_power.o smooth.o
	$(LINK) $(LFLAGS) -lfftw3 $^ -o $@

statistic: statistic.o calc_power.o mean_flux.o smooth.o powerspectrum.o
	$(LINK) $(LFLAGS) -lfftw3 $^ -o $@

%.o: %.c global_vars.h
main.o: main.cpp $(EXT_INC)
read_snapshot.o: read_snapshot.cpp global_vars.h
	$(CXX) $(CFLAGS) -I${GREAD} -c $<

%.o: %.cpp %.h
part_int.o: part_int.cpp part_int.h absorption.h index_table.h

calc_power.o: calc_power.c smooth.o powerspectrum.o

py_module.o: py_module.cpp $(COM_INC)
	$(CXX) $(CFLAGS) -fno-strict-aliasing -DNDEBUG $(PYINC) -c $< -o $@

_spectra_priv.so: py_module.o absorption.o index_table.o part_int.o absorption.o
	$(LINK) $(LFLAGS) -shared $^ -o $@

clean:
	rm -f *.o  extract rescale statistic _spectra_priv.so

dist: Makefile calc_power.c extract_spectra.cpp absorption.cpp py_module.cpp global_vars.h main.cpp mean_flux.c $(COM_INC) powerspectrum.c read_hdf_snapshot.c read_snapshot.cpp rescale.c smooth.c statistic.c statistic.h init.c index_table.cpp
	tar -czf flux_extract.tar.gz $^

btest: test.cpp absorption.o index_table.o
	${LINK} -I${GREAD} ${LFLAGS} -lboost_unit_test_framework $^ -o $@

test: btest
	@./btest

