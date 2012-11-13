#Change this to where you installed GadgetReader
GREAD=${CURDIR}/../GadgetReader
#Python include path
PYINC=-I/usr/local/python/include/python2.6 -I/usr/local/python/x86_64/include/python2.6

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
  CFLAGS +=-O2 -g -c -w1 -openmp -I${GREAD} -fpic
  LINK +=${CXX} -openmp
else
  CFLAGS +=-O2 -g -c -Wall -fopenmp -I${GREAD} -fPIC
  LINK +=${CXX} -openmp $(PRO)
  LFLAGS += -lm -lgomp
endif
#CC = icc -openmp -vec_report0
#CC= gcc -fopenmp -Wall 
OPTS = 
PG = 
OPTS += -DPERIODIC
# Use a periodic spectra
OPTS += -DPECVEL
# Use peculiar velocites 
OPTS += -DVOIGT
# Voigt profiles vs. Gaussian profiles
OPTS += -DHDF5
#Support for loading HDF5 files
OPTS += -DGADGET3 #-DJAMIE
#This is misnamed: in reality it looks for NE instead of NHEP and NHEPP
#OPTS += -DHELIUM
# Enable helium absorption
CFLAGS += $(OPTS) 
CXXFLAGS += $(CFLAGS) -I${GREAD}
COM_INC = parameters.h
#LINK=$(CC)
LFLAGS+=-lfftw3 -lrgad -L${GREAD} -Wl,-rpath,${GREAD} -lhdf5 -lhdf5_hl

.PHONY: all clean dist

all: extract statistic rescale spectra.so

extract: main.o read_snapshot.o read_hdf_snapshot.o extract_spectra.o init.o
	$(LINK) $(LFLAGS) $^ -o $@

rescale: rescale.o powerspectrum.o mean_flux.o calc_power.o smooth.o $(COM_INC)
	$(LINK) $(LFLAGS) $^ -o $@

statistic: statistic.o calc_power.o mean_flux.o smooth.o powerspectrum.o $(COM_INC)
	$(LINK) $(LFLAGS) $^ -o $@

init.o: init.c $(COM_INC)
rescale.o: rescale.c $(COM_INC)
statistic.o: statistic.c $(COM_INC)
read_snapshot.o: read_snapshot.cpp $(COM_INC)
read_hdf_snapshot.o: read_hdf_snapshot.c $(COM_INC)
extract_spectra.o: global_vars.h extract_spectra.c $(COM_INC)
smooth.o:smooth.c
calc_power.o: calc_power.c smooth.o powerspectrum.o 
powerspectrum.o: powerspectrum.c  $(COM_INC)
mean_flux.o: mean_flux.c $(COM_INC)
main.o: main.c global_vars.h $(COM_INC)

py_module.o: py_module.c
	$(CC) $(CFLAGS) -fno-strict-aliasing -DNDEBUG $(PYINC) -c $^ -o $@
spectra.so: py_module.o extract_spectra.o init.o
	$(LINK) -shared $^ -o $@
clean:
	rm -f *.o  extract rescale statistic

dist: Makefile calc_power.c extract_spectra.c global_vars.h main.c mean_flux.c $(COM_INC) powerspectrum.c read_hdf_snapshot.c read_snapshot.cpp rescale.c smooth.c statistic.c statistic.h
	tar -czf flux_extract.tar.gz $^

