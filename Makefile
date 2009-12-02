
CC = gcc
CFLAGS =  -O2 -Wall  -g -pg
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

extract: read_snapshot.o extract_spectra.o readgadget.o Makefile	
	$(CC) $(CFLAGS) -o extract -pg read_snapshot.o -pg extract_spectra.o -pg readgadget.o -lm

read_snapshot.o: read_snapshot.c global_vars.h headers.h parameters.h Makefile
	$(CC) $(CFLAGS) -c read_snapshot.c

readgadget.o: readgadget.c global_vars.h headers.h parameters.h Makefile
	$(CC) $(CFLAGS) -c readgadget.c

extract_spectra.o: extract_spectra.c global_vars.h headers.h parameters.h Makefile
	$(CC) $(CFLAGS) -c extract_spectra.c

clean:
	rm -f *.o  extract
