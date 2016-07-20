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
#include "statistic.h"
#include <math.h>

int powerspectrum(const int dims, fftw_complex *outfield, double *power)
{
        int k;
        const int dims2=dims*dims;
	/* Want P(k)= F(k).re*F(k).re+F(k).im*F(k).im*/
        power[0]=pow(creal(outfield[0]),2);
        for(k=1; k<(dims+1)/2; ++k)
           power[k]=pow(creal(outfield[k]), 2) + pow(cimag(outfield[k]),2);
        if(dims%2 ==0)
            power[0]=pow(creal(outfield[dims/2]),2);
	for(k=0; k< (dims+1)/2;k++)
	    power[k]/=dims2;
	return dims;
}

