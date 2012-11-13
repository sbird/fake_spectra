#include "global_vars.h"
#include <errno.h>
#include <string.h>

int compare_xx(const void *a, const void *b)
{
  if(((sort_los *) a)->priax < (((sort_los *) b)->priax))
    return -1;

  if(((sort_los *) a)->priax > (((sort_los *) b)->priax))
    return +1;

  return 0;
}

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

/*Populate the los table index*/
void populate_sort_los_table(los * los_table, int NumLos, sort_los * sort_los_table, int * nxx)
{
        int i,nother=0;
        /*Make a table with a bit more indirection, so we can sort it*/
        /*Need a pointer to the separate structure for los with iaxis=1*/
        sort_los *sort_los_table_xx;
        for(i=0;i<NumLos;i++){
            if(los_table[i].axis==1){
                  sort_los_table[NumLos-1-*nxx].orig_index=i;
                  sort_los_table[NumLos-1-*nxx].priax=los_table[i].yy;
                  (*nxx)++;
            }else{
                  sort_los_table[nother].orig_index=i;
                  sort_los_table[nother].priax=los_table[i].xx;
                  nother++;
            }
        }
        sort_los_table_xx=sort_los_table+NumLos-*nxx;
        /*Sort the tables: now the table is sorted we can use bsearch to find the element we are looking for*/
        qsort(sort_los_table,NumLos-*nxx,sizeof(sort_los),compare_xx);
        qsort(sort_los_table_xx,*nxx,sizeof(sort_los),compare_xx);

        return;
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
