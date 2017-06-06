#Optional file formats for the C reader.
#Does not affect python.
#Define this to disable HDF5
#HDF5=no
#Define this to disable reading of type 2 snapshots
GREAD=no
#Change this to where you installed GadgetReader
GREADDIR=${CURDIR}/../GadgetReader

CFLAGS +=-O3 -g -Wall -fPIC -ffast-math -fopenmp
LIBS= -lm -lgomp

ifeq ($(HDF5),no)
#Support for loading HDF5 files
CFLAGS += -DNOHDF5
else
LIBS+=-lhdf5 -lhdf5_hl
endif
ifeq ($(GREAD),no)
#Support for loading HDF5 files
CFLAGS += -DNOGREAD
else
IGREAD = -I${GREADDIR}
LIBS+=-lrgad -L${GREADDIR} -Wl,-rpath,${GREADDIR}
endif

COM_INC = fake_spectra/index_table.h fake_spectra/absorption.h fake_spectra/part_int.h

.PHONY: clean test all

all: btest cextract/build cextract/build/extract

btest:
	mkdir btest

cextract/build:
	mkdir cextract/build

cextract/build/extract: cextract/build/main.o cextract/build/read_snapshot.o cextract/build/read_hdf_snapshot.o btest/absorption.o cextract/build/init.o btest/index_table.o btest/part_int.o btest/Faddeeva.o
	$(CXX) $(CFLAGS) $(LIBS) $^ -o $@

cextract/build/statistic: cextract/build/statistic.o cextract/build/calc_power.o cextract/build/mean_flux.o cextract/build/smooth.o cextract/build/powerspectrum.o
	$(CC) $(CFLAGS) -lfftw3 $^ -o $@

cextract/build/%.o: cextract/%.c cextract/global_vars.h cextract/build
	$(CC) $(CFLAGS) ${IGREAD} -c $< -o $@

cextract/build/main.o: cextract/main.cpp cextract/global_vars.h $(COM_INC)
	$(CC) $(CFLAGS) -Ifake_spectra ${IGREAD} -c $< -o $@

cextract/build/read_snapshot.o: cextract/read_snapshot.cpp cextract/global_vars.h
	$(CC) $(CFLAGS) ${IGREAD} -c $< -o $@

btest/%.o: fake_spectra/%.cpp fake_spectra/%.h
	$(CXX) $(CFLAGS) -c $< -o $@

btest/absorption.o: fake_spectra/absorption.cpp fake_spectra/absorption.h fake_spectra/singleabs.h
btest/part_int.o: fake_spectra/part_int.cpp fake_spectra/part_int.h fake_spectra/absorption.h fake_spectra/index_table.h

clean:
	rm -f cextract/build/* btest/*

btest/btest: fake_spectra/test.cpp btest/absorption.o btest/index_table.o btest/Faddeeva.o
	$(CXX) $(CFLAGS) $^ -lboost_unit_test_framework -o $@

test: btest btest/btest
	@./btest/btest

