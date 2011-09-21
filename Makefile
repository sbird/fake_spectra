#Change this to where you installed GadgetReader
GREAD=${CURDIR}/../GadgetReader

#CC = icc -openmp -vec_report0
#CC= gcc -fopenmp -Wall 
CFLAGS =  -O2  -g -fopenmp -Wall
OPTS = 
PG = 
OPTS += -DPERIODIC
# Use a periodic spectra
OPTS += -DPECVEL
# Use peculiar velocites 
OPTS += -DVOIGT
# Voigt profiles vs. Gaussian profiles
#OPTS += -DGADGET3
#Gadget III has slightly different block headers
#OPTS += -DHELIUM
# Enable helium absorption
CFLAGS += $(OPTS) 
CXXFLAGS += $(CFLAGS) -I${GREAD}
COM_INC = parameters.h Makefile
FFTW =-lfftw3
#LINK=$(CC)
LINK=$(CXX) -lm -lgomp -lfftw3 -lrgad -L$(FFTW) -L${GREAD} -Wl-rpath,${GREAD}

.PHONY: all clean

all: extract statistic

extract: main.o read_snapshot.o extract_spectra.o Makefile
	$(LINK) $(CFLAGS) -o extract $(PG) main.o $(PG) read_snapshot.o $(PG) extract_spectra.o $(PG) 

rescale: rescale.c powerspectrum.o mean_flux.o calc_power.o smooth.o $(COM_INC)
	$(CC) $(CFLAGS) rescale.c $(FFTW) powerspectrum.o mean_flux.o calc_power.o smooth.o -o rescale 

statistic: statistic.c calc_power.o mean_flux.o smooth.o $(COM_INC)
	$(CC) $(CFLAGS) statistic.c  powerspectrum.o mean_flux.o calc_power.o smooth.o -o statistic $(FFTW)

read_snapshot.o: read_snapshot.cpp $(COM_INC)
extract_spectra.o: global_vars.h extract_spectra.c $(COM_INC)
smooth.o:smooth.c
calc_power.o: calc_power.c smooth.o powerspectrum.o 
powerspectrum.o: powerspectrum.c  $(COM_INC)
mean_flux.o: mean_flux.c $(COM_INC)
main.o: main.c global_vars.h $(COM_INC)

clean:
	rm -f *.o  extract


