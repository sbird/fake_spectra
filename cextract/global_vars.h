#ifndef GLOBAL_VARS_H
#define GLOBAL_VARS_H

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

/* Spectrum data set size     */

struct particle_data
{
  float *Pos;
  float *Vel;
  float *Mass;
  float *U, *Ne;
  float *temp;
  float *h;
  float *fraction;
};
typedef struct particle_data pdata;

/*Functions to allocate memory.*/
void help(void);

#ifndef NOGREAD
/*Functions to load Gadget-2 snapshots.*/
int64_t load_snapshot(const char *fname, int64_t StartPart,pdata* P);
int load_header(const char *fname,double  *atime, double *redshift, double * Hz, double *box100, double *h100);
#endif

#ifdef __cplusplus
extern "C"{
#endif
int alloc_parts(pdata* P, int np);
void free_parts(pdata* P);

#ifndef NOHDF5
int load_hdf5_header(const char *infname, double  *atime, double *redshift, double * Hz, double *box100, double *h100);
int load_hdf5_snapshot(const char *ffname, pdata *P, int fileno);
#endif
void populate_los_table(double * cofm, int * axis, int NumLos, char * ext_table, double box);
#ifdef __cplusplus
}
#endif

#endif
