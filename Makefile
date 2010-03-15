
#CC = icc -openmp -vec_report0
CC= gcc -fopenmp -Wall 
CFLAGS =  -O2  -g
OPTS = 
PG = 
OPTS += -DPERIODIC
# Use a periodic spectra
OPTS += -DPECVEL
# Use peculiar velocites 
OPTS += -DVOIGT 
# Voigt profiles vs. Gaussian profiles
#OPTS += -DHELIUM
# Enable helium absorption
OPTS += -DNO_RESCALE_FLUX
#Don't rescale to the observed mean flux
CFLAGS += $(OPTS)
COM_INC = global_vars.h parameters.h Makefile
LINK=$(CC) -lm -lsrfftw -lsfftw -lgomp -L/data/store/spb41/apps/fftw/lib

.PHONY: all clean

all: extract

extract-icc: main.o read_snapshot.o extract_spectra.o readgadget.o powerspectrum.o mean_flux.o Makefile	
	$(LINK) $(CFLAGS) -o extract-icc $(PG) main.o $(PG) read_snapshot.o $(PG) extract_spectra.o $(PG) readgadget.o $(PG) powerspectrum.o $(PG) mean_flux.o

extract: main.o read_snapshot.o extract_spectra.o readgadget.o powerspectrum.o mean_flux.o Makefile	
	$(LINK) $(CFLAGS) -o extract $(PG) main.o $(PG) read_snapshot.o $(PG) extract_spectra.o $(PG) readgadget.o $(PG) powerspectrum.o $(PG) mean_flux.o
read_snapshot.o: read_snapshot.c $(COM_INC)
	$(CC) $(CFLAGS) -c read_snapshot.c

readgadget.o: readgadget.c $(COM_INC)
	$(CC) $(CFLAGS) -c readgadget.c

extract_spectra.o: extract_spectra.c $(COM_INC)
	$(CC) $(CFLAGS) -c extract_spectra.c

powerspectrum.o: powerspectrum.c  $(COM_INC)
	$(CC) $(CFLAGS) -c powerspectrum.c

mean_flux.o: mean_flux.c $(COM_INC)
	$(CC) $(CFLAGS) -c mean_flux.c

main.o: main.c $(COM_INC)
	$(CC) $(CFLAGS) -c main.c

clean:
	rm -f *.o  extract
