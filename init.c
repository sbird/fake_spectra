#include "global_vars.h"
#include <errno.h>
#include <string.h>

/* Populate the line of sight table, either by random numbers or with some external input. */
void populate_los_table(los * los_table, int NumLos, char * ext_table, double box)
{
        FILE * fh;
        int lines=0;
        int axis;
        float xx, yy, zz;
        /*If we have a file path, load the sightline table from there*/
        if(ext_table){
                const double boxm=box/1000.0;
                if(!(fh=fopen(ext_table, "r")))
                {
                        fprintf(stderr, "Error opening %s: %s\n",ext_table, strerror(errno));
                        exit(3);
                }
                while(lines < NumLos){
                        if(EOF == fscanf(fh, "%d %f %f %f\n", &axis, &xx, &yy, &zz)){
                                fprintf(stderr, "Error reading table: %s. Possibly file is truncated?\n",strerror(errno));
                                exit(3);
                        }
                        if(axis > 3 || axis <0){
                                fprintf(stderr, "Line %d of gives axis %d, which is silly.\n", lines+1, axis);
                                exit(3);
                        }
                        if (xx > boxm || xx < 0 ||
                           yy > boxm || yy < 0 || zz > boxm || zz <0 ){
                                fprintf(stderr, "Line %d of LOS table is: %d %f %f %f, which is silly for boxm %f.\n", lines+1, axis, xx, yy, zz, boxm);
                                exit(3);
                        }
                        los_table[lines].axis=axis;
                        los_table[lines].xx=xx*1000;
                        los_table[lines].yy=yy*1000;
                        los_table[lines].zz=zz*1000;
                        lines++;
                }
        }
        else{
                srand48(241008); /* random seed generator */
                for(lines=0; lines<NumLos; lines++)
                {
                        do	
                        	axis = (int)(drand48()*4);
                        while (axis == 0 || axis==4); 
                        los_table[lines].axis=axis;
                        los_table[lines].xx=drand48()*box;
                        los_table[lines].yy=drand48()*box;
                        los_table[lines].zz=drand48()*box;
                }
        }
}

/*****************************************************************************/
int InitLOSMemory(interp* species,int NumLos, int nbins)
{  
  (*species).rho        = (double *) calloc((NumLos * nbins) , sizeof(double));
  (*species).veloc        = (double *) calloc((NumLos * nbins) , sizeof(double));
  (*species).temp   = (double *) calloc((NumLos * nbins) , sizeof(double));
  if(!(*species).rho || !(*species).veloc || !(*species).temp)
      return 1;
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
void FreeLOSMemory(interp * species)
{  
  free((*species).rho);
  free((*species).veloc);
  free((*species).temp);
}

/*Note this assumes only one species*/
int alloc_parts(pdata* P, int np)
{
    return ((*P).Vel=(float *)malloc(np*3*sizeof(float))) &&
    ((*P).Pos=(float *)malloc(np*3*sizeof(float))) &&
     ((*P).Mass=(float *) malloc(np*sizeof(float))) &&
    ((*P).U=(float *)malloc(np*sizeof(float))) &&
    ((*P).fraction=(float *)malloc(np*sizeof(float))) &&
    ((*P).Ne=(float *)malloc(np*sizeof(float))) &&
    ((*P).h=(float *)malloc(np*sizeof(float)));
}

void free_parts(pdata* P)
{
    free((*P).Vel);
    free((*P).Pos);
    free((*P).Mass);
    free((*P).U);
    free((*P).fraction);
    free((*P).Ne);
    free((*P).h);
    return;
}
