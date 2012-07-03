#ifndef GLOBAL_VARS_H
#define GLOBAL_VARS_H

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>


struct particle_data 
{
  float *Pos;
  float *Vel;
  float *Mass;
  float *U, *NH0, *Ne, *h;
#ifdef HELIUM
  float *NHep;
#endif
};
typedef struct particle_data pdata;

struct _interp
{
   double *rho;
   double *temp;
   double *veloc;
};
typedef struct _interp interp;

/*Are we arepo?*/
static int arepo;

/*Allocate and free memory for the particle tables*/
#ifdef __cplusplus
extern "C"
#endif
int alloc_parts(pdata* P, int np);
#ifdef __cplusplus
extern "C"
#endif
void free_parts(pdata* P);

/*Structure for storing a sightline*/
struct _los
{
        int axis;
        float xx;
        float yy;
        float zz;
};
typedef struct _los los;

/*Structure and comparison function for sorted list of los*/
typedef struct _sort_los
{
        int orig_index;
        double priax;
        /*This is xx, unless iaxis=1, in which case it is yy*/
} sort_los;

/*Pointers to arrays to use in SPH_interpolation*/
#ifndef RAW_SPECTRA
#include "statistic.h"
#endif

/*Functions to allocate memory.*/
int InitLOSMemory(interp * species, int NumLos);
void FreeLOSMemory(interp * species);
int WriteLOSData(interp* species,double * tau, int NumLos,FILE * output);

void help(void);

#ifdef __cplusplus
extern "C"{
#endif
int64_t load_snapshot(char *fname, int64_t StartPart,int64_t MaxRead,pdata* P, double *omegab);
int load_header(char *fname,double  *atime, double *redshift, double * Hz, double *box100, double *h100);
#ifdef __cplusplus
}
#endif
#ifdef HDF5
int load_hdf5_header(char *infname, double  *atime, double *redshift, double * Hz, double *box100, double *h100);
int load_hdf5_snapshot(char *ffname, pdata *P, double *omegab, int fileno);
int find_first_hdf_file(const char *infname, char *fname);
#endif
void populate_los_table(los * los_table, int NumLos, sort_los * sort_los_table, int * nxx, char * ext_table, double box);

#ifndef HELIUM
void SPH_Interpolation(double * rhoker_H, interp * H1, const int Particles, const int NumLos,const double boxsize, const los *los_table,const sort_los *sort_los_table,const int nxx, const pdata *P);
void Compute_Absorption(double * tau_H1, double *rhoker_H,interp * H1, const double Hz, const double h100, const double box100, const double atime, const double omegab);
#else
void Compute_Absorption(double * tau_H1, double *rhoker_H, interp * H1,double * tau_He2,interp * He2, const double Hz, const double h100, const double box100, const double atime, const double omegab);
void SPH_Interpolation(double * rhoker_H, interp * H1, interp * He2, const int Particles, const int NumLos,const double boxsize, const los *los_table,const sort_los *sort_los_table,const int nxx, const pdata *P);
#endif

#endif
