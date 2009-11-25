
CC = icc
CFLAGS =  -O -w  -i64 -ip 

extract: read_snapshot.o extract_spectra.o	
	$(CC) $(CFLAGS) -o extract read_snapshot.o extract_spectra.o -lm

read_snapshot.o: read_snapshot.c global_vars.h headers.h parameters.h
	$(CC) $(CFLAGS) -c read_snapshot.c

extract_spectra.o: extract_spectra.c global_vars.h headers.h parameters.h
	$(CC) $(CFLAGS) -c extract_spectra.c

clean:
	rm -f *.o  extract
