#ifndef GLOBAL_VARS_H
#define GLOBAL_VARS_H

#include <stdlib.h>
#include <stdio.h>

/*Gadget particle header*/
struct io_header_1
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  int      flag_stellarage;
  int      flag_metals;
  char     fill[88];  /* fills to 256 Bytes */
};
typedef struct io_header_1 gadget_header;

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

/*Allocate and free memory for the particle tables*/
int alloc_parts(pdata* P, int np);
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

/*Pointers to arrays to use in SPH_interpolation*/
#ifndef RAW_SPECTRA
#include "statistic.h"
#endif

/*Functions to allocate memory.*/
int InitLOSMemory(interp * species, int NumLos);
void FreeLOSMemory(interp * species);
int WriteLOSData(interp* species,double * tau, int NumLos,FILE * output);

#define int_blk int64_t
void swap_Nbyte(char *data,int n,int m);
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream); 
int_blk find_block(FILE *fd,char *label);
int_blk read_gadget_float(float *data,char *label,int offset, int read,FILE *fd, int old);
/* The final argument, if one, means it will attempt to read an old format file
 * It may not work, due to the wide variation in GADGET type one files. */
int_blk read_gadget_float3(float *data,char *label,int offset, int read, FILE *fd, int old);
int read_gadget_head(gadget_header *out_header, FILE *fd, int old);
void help(void);

int load_snapshot(char *fname, int files, int old, pdata* P,
  double  *atime, double *redshift, double * Hz, double *box100, double *h100, double *omegab);
void populate_los_table(los *los_table, int NumLos, char *ext_table, double box);

#ifndef HELIUM
void SPH_Interpolation(double * rhoker_H, interp * H1, const int Particles, const int NumLos,const double boxsize, const los *los_table, const pdata *P);
void Compute_Absorption(double * tau_H1, double *rhoker_H,interp * H1, const double Hz, const double h100, const double box100, const double atime, const double omegab);
#else
void Compute_Absorption(double * tau_H1, double *rhoker_H, interp * H1,double * tau_He2,interp * He2, const double Hz, const double h100, const double box100, const double atime, const double omegab);
void SPH_Interpolation(double * rhoker_H, interp * H1, interp * He2, const int Particles, const int NumLos,const double boxsize, const los *los_table, const pdata *P);
#endif

#endif
