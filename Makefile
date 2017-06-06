#Optional file formats for the C reader.
#Does not affect python.
#Define this to disable HDF5
#HDF5=no
#Define this to disable reading of type 2 snapshots
GREAD=no
#Change this to where you installed GadgetReader
GREADDIR=${CURDIR}/../GadgetReader

# Voigt profiles vs. Gaussian profiles
CFLAGS += -DVOIGT -Ifake_spectra

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

COM_INC = fake_spectra/index_table.h fake_spectra/absorption.h fake_spectra/part_int.h
CXXFLAGS += $(CFLAGS)

.PHONY: clean test all

all: btest cextract/build cextract/build/extract

btest:
	mkdir btest

cextract/build:
	mkdir cextract/build

cextract/build/extract: cextract/build/main.o cextract/build/read_snapshot.o cextract/build/read_hdf_snapshot.o btest/absorption.o cextract/build/init.o btest/index_table.o btest/part_int.o btest/Faddeeva.o
	$(LINK) $(LFLAGS) $(LIBS) $^ -o $@

cextract/build/statistic: cextract/build/statistic.o cextract/build/calc_power.o cextract/build/mean_flux.o cextract/build/smooth.o cextract/build/powerspectrum.o
	$(LINK) $(LFLAGS) -lfftw3 $^ -o $@

cextract/build/%.o: cextract/%.c cextract/global_vars.h cextract/build
	$(CC) $(CFLAGS) ${IGREAD} -c $< -o $@

cextract/build/main.o: cextract/main.cpp cextract/global_vars.h $(COM_INC)
	$(CC) $(CFLAGS) ${IGREAD} -c $< -o $@

cextract/build/read_snapshot.o: cextract/read_snapshot.cpp cextract/global_vars.h
	$(CC) $(CFLAGS) ${IGREAD} -c $< -o $@

btest/%.o: fake_spectra/%.cpp fake_spectra/%.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

btest/absorption.o: fake_spectra/absorption.cpp fake_spectra/absorption.h fake_spectra/singleabs.h
btest/part_int.o: fake_spectra/part_int.cpp fake_spectra/part_int.h fake_spectra/absorption.h fake_spectra/index_table.h

clean:
	rm -f cextract/build/* btest/*

btest/btest: fake_spectra/test.cpp btest/absorption.o btest/index_table.o btest/Faddeeva.o
	${LINK} ${LFLAGS} $^ -lboost_unit_test_framework -o $@

test: btest btest/btest
	@./btest/btest

