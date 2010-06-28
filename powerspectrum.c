/* Copyright (c) 2009, Simeon Bird <spb41@cam.ac.uk>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE. */
#include <drfftw.h>
#include "global_vars.h"
#include <omp.h>

/* Computes the one-dimensional normalized power spectrum*/
/*Little macro to work the storage order of the FFT.*/
#define KVAL(n) ((n)<=dims/2 ? (n) : ((n)-dims))

int powerspectrum(const int dims, double *field, double *power)
{
        fftw_real *outfield;
        int k;
        const int dims2=dims*dims;
	if(sizeof(fftw_real) != sizeof(double))
	{
		fprintf(stderr, "sizeof fftw_real:%zd fftw_complex: %zd, double: %zd\n",sizeof(fftw_real), sizeof(fftw_complex), sizeof(double));
		fprintf(stderr, "fftw_real is not a double. Perhaps you linked the wrong library?\n");
		exit(1);
	}
	outfield=malloc(dims*sizeof(fftw_real));
	if(!outfield){
			  fprintf(stderr, "Error allocating memory for outfield!\n");
			  exit(1);
	}
	rfftw_one(pl, &field[0], outfield);
	/* Want P(k)= F(k).re*F(k).re+F(k).im*F(k).im
	 * FFTW has a strange format for output array:
    * F(0).re F(1).re ..F(n/2).re F((n+1)/2-1).im ...F(0).im */
        power[0]=outfield[0]*outfield[0];
        for(k=1; k<(dims+1)/2; ++k)
           power[k]=outfield[k]*outfield[k]+outfield[dims-k]*outfield[dims-k];
        if(dims%2 ==0)
            power[dims/2]=outfield[dims/2]*outfield[dims/2];
	for(k=0; k< (dims+1)/2;k++)
	    power[k]/=dims2;
        free(outfield);
	return dims;
}

