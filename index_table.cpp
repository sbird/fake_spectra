#include "global_vars.h"
#include <math.h>

int compare_xx(const void *a, const void *b)
{
  if(((sort_los *) a)->priax < (((sort_los *) b)->priax))
    return -1;

  if(((sort_los *) a)->priax > (((sort_los *) b)->priax))
    return +1;

  return 0;
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
        qsort(sort_los_table,nother,sizeof(sort_los),compare_xx);
        qsort(sort_los_table_xx,*nxx,sizeof(sort_los),compare_xx);

        return;
}


/*This implements binary search for a given xx.
 * Returns the index of the last element where xx >= priax */
int find_index(double xx, const sort_los* sort_los_table, const int NumLos)
{
        int low,high,mid;
        low=0;
        high=NumLos-1;
        if(xx < sort_los_table[0].priax)
                return 0;
        while(high - low > 1)
        {
            mid = (high + low) / 2;
            if(xx < sort_los_table[mid].priax)
                high = mid;
            else
                low = mid;
        }
        return high;
}

int get_near_lines_2nd_axis(const double xx,const double yy,const double zz,const double h4, const double boxsize,const sort_los *sort_los_table, const los *los_table, int *index_nr_lines, double *dr2_lines, const int low, const int high)
{
      int ind,num_nr_lines=0;
      for(ind=low;ind<=high;ind++)
      {
          double dr,dr2;
          const int iproc=sort_los_table[ind].orig_index;
          /*Load a sightline from the table.*/
          const int iaxis = los_table[iproc].axis;
          double xproj,yproj,zproj;

          xproj = los_table[iproc].xx;
          yproj = los_table[iproc].yy;
          zproj = los_table[iproc].zz;

          /*    Distance to projection axis */
          if (iaxis == 1)
            dr = fabs(yy-yproj);
          else
            dr = fabs(xx-xproj);

          if(dr > 0.5*boxsize)
                  dr = boxsize - dr; /* Keep dr between 0 and box/2 */

          dr2 = dr*dr;

          if (iaxis == 3)
            dr = fabs(yy - yproj);
          else
            dr = fabs(zz - zproj);

          if (dr > 0.5*boxsize)
            dr = boxsize - dr; /* between 0 and box/2 */

          dr2 = dr2 + (dr*dr);

          /*If close in the second coord, save line*/
          if (dr2 <= h4){
                  index_nr_lines[num_nr_lines]=iproc;
                  dr2_lines[num_nr_lines]=dr2;
                  num_nr_lines++;
          }
      }
      return num_nr_lines;
}

/*This function takes a particle position and returns a list of the indices of lines near it in index_nr_lines
 * Near is defined as: dx^2+dy^2 < 4h^2 */
int get_list_of_near_lines(const double xx,const double yy,const double zz,const double hh, const double boxsize,const los *los_table, const int NumLos,const sort_los* sort_los_table,const int nxx, int *index_nr_lines, double *dr2_lines)
{
      const double h4 = 4.*hh*hh;           /* 2 smoothing lengths squared */
      int low,high;
      int num_nr_lines=0;
      double ff;
      /*Need a pointer to the separate structure for los with iaxis=1*/
      const sort_los *sort_los_table_xx;
      sort_los_table_xx=sort_los_table+NumLos-nxx;
      if(nxx < NumLos){
        /*Now find the elements where dr < 2 hh, wrapping with respect to boxsize*/
        /* First find highest index where xx + 2 hh > priax */
        ff=xx+2*hh;
        if(ff > boxsize)
                ff-=boxsize;
        high=find_index(ff,sort_los_table,NumLos-nxx);
        /* Now find lowest index in what remains where xx - 2 hh < priax */
        ff=xx-2*hh;
        if(ff < 0)
                ff+=boxsize;
        low=find_index(ff,sort_los_table,NumLos-nxx);
        /*This should be the case unless wrapping has occurred*/
        if(low <= high)
          num_nr_lines+=get_near_lines_2nd_axis(xx,yy,zz,h4, boxsize,sort_los_table, los_table, index_nr_lines+num_nr_lines, dr2_lines, low, high);
        else{
          num_nr_lines+=get_near_lines_2nd_axis(xx,yy,zz,h4, boxsize,sort_los_table, los_table, index_nr_lines+num_nr_lines, dr2_lines, 0, high);
          num_nr_lines+=get_near_lines_2nd_axis(xx,yy,zz,h4, boxsize,sort_los_table, los_table, index_nr_lines+num_nr_lines, dr2_lines, low, NumLos-nxx);
        }
      }
      if(nxx > 0){
        /*Do the same thing with the table where iaxis=1*/
        /*Now find the elements where dr < 2 hh, wrapping with respect to boxsize*/
        /* First find highest index where xx + 2 hh >= priax */
        ff=yy+2*hh;
        if(ff > boxsize)
                ff-=boxsize;
        high=find_index(ff,sort_los_table_xx,nxx);
        /* Now find highest index in what remains where xx - 2 hh >= priax */
        ff=yy-2*hh;
        if(ff < 0)
                ff+=boxsize;
        low=find_index(ff,sort_los_table_xx,nxx);
        /*This should be the case unless wrapping has occurred*/
        if(low <= high)
          num_nr_lines+=get_near_lines_2nd_axis(xx,yy,zz,h4, boxsize,sort_los_table_xx, los_table, index_nr_lines+num_nr_lines, dr2_lines, low, high);
        else{
          num_nr_lines+=get_near_lines_2nd_axis(xx,yy,zz,h4, boxsize,sort_los_table_xx, los_table, index_nr_lines+num_nr_lines, dr2_lines, 0, high);
          num_nr_lines+=get_near_lines_2nd_axis(xx,yy,zz,h4, boxsize,sort_los_table_xx, los_table, index_nr_lines+num_nr_lines, dr2_lines, low, nxx);
        }
      }
      return num_nr_lines;
}
