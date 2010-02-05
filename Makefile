
CC = gcc
CFLAGS =  -O2 -Wall  -g -fopenmp
OPTS = 
OPTS += -DPERIODIC
# Use a periodic spectra
OPTS += -DPECVEL
# Use peculiar velocites 
OPTS += -DVOIGT 
# Voigt profiles vs. Gaussian profiles
#OPTS += -DHELIUM
# Enable helium absorption
CFLAGS += $(OPTS)
LINK=gcc -lm -lsrfftw -lsfftw -lgomp -L/data/store/spb41/apps/fftw/lib

extract: read_snapshot.o extract_spectra.o readgadget.o powerspectrum.o mean_flux.o Makefile	
	$(LINK) $(CFLAGS) -o extract read_snapshot.o extract_spectra.o readgadget.o powerspectrum.o mean_flux.o

read_snapshot.o: read_snapshot.c global_vars.h headers.h parameters.h Makefile
	$(CC) $(CFLAGS) -c read_snapshot.c

readgadget.o: readgadget.c global_vars.h headers.h parameters.h Makefile
	$(CC) $(CFLAGS) -c readgadget.c

extract_spectra.o: extract_spectra.c global_vars.h headers.h parameters.h Makefile
	$(CC) $(CFLAGS) -c extract_spectra.c

powerspectrum.o: powerspectrum.c
	$(CC) $(CFLAGS) -c powerspectrum.c

mean_flux.o: mean_flux.c
	$(CC) $(CFLAGS) -c mean_flux.c
clean:
	rm -f *.o  extract
