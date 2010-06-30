
CC = icc -openmp -vec_report0
#CC= gcc -fopenmp -Wall 
CFLAGS =  -O2  -g
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
OPTS += -DRAW_SPECTRA
# Output data file with optical depths rather than flux power spectrum.
#OPTS += -DHELIUM
# Enable helium absorption
CFLAGS += $(OPTS)
COM_INC = parameters.h Makefile
FFTW =-ldrfftw -ldfftw
LINK=$(CC)
#LINK=$(CC) -lm -lgomp -lsrfftw -lsfftw  -L$(FFTW)

.PHONY: all clean

all: extract

extract: main.o read_snapshot.o extract_spectra.o readgadget.o Makefile
	# powerspectrum.o mean_flux.o calc_power.o smooth.o
	$(LINK) $(CFLAGS) -o extract $(PG) main.o $(PG) read_snapshot.o $(PG) extract_spectra.o $(PG) readgadget.o 
	#$(PG) powerspectrum.o $(PG) mean_flux.o calc_power.o smooth.o $(FFTW)

rescale: rescale.c powerspectrum.o mean_flux.o calc_power.o smooth.o $(COM_INC)
	$(CC) $(CFLAGS) rescale.c $(FFTW) powerspectrum.o mean_flux.o calc_power.o smooth.o -o rescale 

statistic: statistic.c calc_power.o $(COM_INC)
	$(CC) $(CFLAGS) statistic.c  powerspectrum.o mean_flux.o calc_power.o smooth.o -o statistic $(FFTW)

read_snapshot.o: read_snapshot.c $(COM_INC)
readgadget.o: readgadget.c $(COM_INC)
extract_spectra.o: global_vars.h extract_spectra.c $(COM_INC)
smooth.o:smooth.c
calc_power.o: calc_power.c smooth.o powerspectrum.o 
powerspectrum.o: powerspectrum.c  $(COM_INC)
mean_flux.o: mean_flux.c $(COM_INC)
main.o: main.c global_vars.h $(COM_INC)

clean:
	rm -f *.o  extract


