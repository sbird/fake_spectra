#ifndef GLOBAL_VARS_H
#define GLOBAL_VARS_H

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include "types.h"

/*Functions to allocate memory.*/
void help(void);


#ifdef __cplusplus
extern "C"{
#endif
int InitLOSMemory(interp * species, int NumLos, int nbins);
void FreeLOSMemory(interp * species);
int WriteLOSData(interp* species,double * tau, int NumLos,FILE * output);

/*Allocate and free memory for the particle tables*/
int alloc_parts(pdata* P, int np);
void free_parts(pdata* P);
int64_t load_snapshot(const char *fname, int64_t StartPart,int64_t MaxRead,pdata* P, double *omegab);
int load_header(const char *fname,double  *atime, double *redshift, double * Hz, double *box100, double *h100);
#ifdef HDF5
int load_hdf5_header(const char *infname, double  *atime, double *redshift, double * Hz, double *box100, double *h100);
int load_hdf5_snapshot(const char *ffname, pdata *P, double *omegab, int fileno);
#endif
void populate_los_table(los * los_table, int NumLos, char * ext_table, double box);
#ifdef __cplusplus
}
#endif

#endif
