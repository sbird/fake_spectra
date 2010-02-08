/* Copyright (c) 2009, Simeon Bird <spb41@cam.ac.uk>
 *               Based on code (c) 2005 by J. Bolton
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global_vars.h"
#include "parameters.h"

/* Here we load a snapshot file. It can be distributed
 * onto several files (for files>1) */
/**********************************************************************/
int main(int argc, char **argv)
{
  int Npart, NumLos, files;
  if(argc<4)
  {
    printf("Usage: ./extract NUMLOS NUMFILES base_filename\n");
    exit(99);
  }
  NumLos=atoi(argv[1]);
  files=atoi(argv[2]);


  if(NumLos <=0 || files <=0)
  {
          printf("Need NUMLOS >0\n");
          exit(99);
  }
  pl=rfftw_create_plan(NBINS,FFTW_REAL_TO_COMPLEX, FFTW_MEASURE | FFTW_THREADSAFE);
  Npart=load_snapshot(argv[3], files);
  InitLOSMemory(NumLos);
  if(!PARTTYPE)
    SPH_interpolation(NumLos,Npart);
  free(P);
  FreeLOSMemory();
  fftw_destroy_plan(pl);
  return 0;
}
/**********************************************************************/
