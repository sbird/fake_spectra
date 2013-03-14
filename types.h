#ifndef TYPES_H
#define TYPES_H

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

/*Structure for storing a sightline*/
struct _los
{
        int axis;
        float xx;
        float yy;
        float zz;
};
typedef struct _los los;

#endif
