#ifndef GLOBAL_VARS_H
#define GLOBAL_VARS_H

/*Plan needs to be generated only once,
 * and used by each thread.*/
#include <srfftw.h>
fftw_plan pl;

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

double  atime, redshift, omega0, omegaL, box100, h100, omegab;

/*Pointers to arrays to use in SPH_interpolation*/
#ifdef RAW_SPECTRA
double *Delta,*n_H1,*veloc_H1,*temp_H1;
#endif
double *tau_H1, *posaxis,*velaxis;
float *flux_power;
#ifdef HELIUM
double *n_He2,*veloc_He2,*temp_He2,*tau_He2;
#endif

/*Functions to allocate memory.*/
void InitLOSMemory(int NumLos);
void FreeLOSMemory(void);

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
/* These functions do the work*/
int powerspectrum(const int dims, float *field, float *power);
double mean_flux(double * tau, double nbins, double obs_flux, double tol);
void calc_power_spectra(float *flux_power, double *tau_H1,double scale, int NumLos);

int load_snapshot(char *fname, int files, int old, pdata* P);
void SPH_interpolation(int NumLos, int Ntype, los *los_table, pdata* P);
void populate_los_table(los *los_table, int NumLos, char *ext_table, double box);

#endif
