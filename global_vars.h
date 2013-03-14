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
  float *U, *Ne, *h;
  float *fraction;
};
typedef struct particle_data pdata;

struct _interp
{
   double *rho;
   double *temp;
   double *veloc;
};
typedef struct _interp interp;

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

/*Functions to allocate memory.*/
int InitLOSMemory(interp * species, int NumLos, int nbins);
void FreeLOSMemory(interp * species);
int WriteLOSData(interp* species,double * tau, int NumLos,FILE * output);

void help(void);

#ifdef __cplusplus
extern "C"{
#endif
int64_t load_snapshot(char *fname, int64_t StartPart,int64_t MaxRead,pdata* P, double *omegab);
int load_header(char *fname,double  *atime, double *redshift, double * Hz, double *box100, double *h100);

void populate_sort_los_table(los * los_table, int NumLos, sort_los * sort_los_table, int * nxx);
int get_list_of_near_lines(const double xx,const double yy,const double zz,const double hh, const double boxsize,const los *los_table, const int NumLos,const sort_los* sort_los_table,const int nxx, int *index_nr_lines, double *dr2_lines);

#ifdef __cplusplus
}
#endif
#ifdef HDF5
int load_hdf5_header(char *infname, double  *atime, double *redshift, double * Hz, double *box100, double *h100);
int load_hdf5_snapshot(char *ffname, pdata *P, double *omegab, int fileno);
int find_first_hdf_file(const char *infname, char *fname);
#endif
void populate_los_table(los * los_table, int NumLos, char * ext_table, double box);
void Compute_Absorption(double * tau_H1, double * rho, double * veloc, double * temp, const int nbins, const double Hz, const double h100, const double box100, const double atime, const double lambda_lya, const double gamma_lya, const double fosc_lya, const double mass);

void SPH_Interpolation(double * rhoker_H, interp * species, const int nspecies, const int nbins, const int Particles, const int NumLos,const double boxsize, const los *los_table,const sort_los *sort_los_table,const int nxx, const pdata *P);
void Rescale_Units(double * rho, double * veloc, double * temp, const int nbins, const double h100, const double atime);
void Convert_Density(double * rhoker_H, double * rho, const double h100, const double atime, const double omegab);

#endif
