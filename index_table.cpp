#include "global_vars.h"
#include "index_table.h"
#include <map>
#include <math.h>
#include <iostream>

// Construct the index tables as maps, which are automatically sorted.
IndexTable::IndexTable(const los * los_table_i, const int NumLos_i, const double box):
  los_table(los_table_i), NumLos(NumLos_i), boxsize(box)
{
        for(int i=0;i<NumLos;i++){
            if(los_table[i].axis==1)
                index_table_xx.insert(std::pair<const double, const int>(los_table[i].yy, i));
            else
                index_table.insert(std::pair<const double, const int>(los_table[i].xx, i));
        }
        if(index_table_xx.size() + index_table.size() != (unsigned int) NumLos){
            std::cerr << "Did not add all elements. xx index: "<<index_table_xx.size()<<" other index: "<<index_table.size()<<std::endl;
        }
        return;
}

/*Returns a std::map of lines close to the coordinates xx, yy, zz.
 * the key is the line index, and the value is the distance from the two axes not projected along*/
void IndexTable::get_nearby_from_range(std::multimap<const double, const int>::iterator low, std::multimap<const double, const int>::iterator high, std::map<int, double>& nearby, const float pos[], const double hh)
{
      for(std::multimap<const double, const int>::iterator it = low; it != high;++it)
      {
          const int iproc=it->second;
          /*If close in the second coord, save line*/

          double dr2 = calc_dr2(iproc, pos);
          if (dr2 <= 4*hh*hh){
                  nearby[iproc]=dr2;
          }
      }
}

double IndexTable::calc_dr2(const int iproc, const float pos[])
{
    double dr, dr2;
    /*Load a sightline from the table.*/
    const int iaxis = los_table[iproc].axis;
    double xproj,yproj,zproj;

    xproj = los_table[iproc].xx;
    yproj = los_table[iproc].yy;
    zproj = los_table[iproc].zz;

    /*    Distance to projection axis */
    if (iaxis == 1)
      dr = fabs(pos[1]-yproj);
    else
      dr = fabs(pos[0]-xproj);

    if(dr > 0.5*boxsize)
            dr = boxsize - dr; /* Keep dr between 0 and box/2 */

    dr2 = dr*dr;

    if (iaxis == 3)
      dr = fabs(pos[1] - yproj);
    else
      dr = fabs(pos[2] - zproj);

    if (dr > 0.5*boxsize)
      dr = boxsize - dr; /* between 0 and box/2 */

    dr2 = dr2 + (dr*dr);
    return dr2;
}

void IndexTable::get_nearby(float first, std::multimap<const double, const int>& sort_los, std::map<int, double>& nearby, const float pos[], const double hh)
{
      /*Now find the elements where dr < 2 hh, wrapping with respect to boxsize*/
      /* First find highest index where xx + 2 hh > priax */
      double ffp=first+2*hh;
      if(ffp > boxsize)
         ffp-=boxsize;
      /* Now find lowest index in what remains where xx - 2 hh < priax */
      double ffm=first-2*hh;
      if(ffm < 0)
         ffm+=boxsize;
      std::multimap<const double, const int>::iterator low,high;
      //An iterator to the first element not less than ffm
      low=sort_los.lower_bound(ffm);
      //An iterator to the first element greater than ffp
      high=sort_los.lower_bound(ffp);
      //If periodic wrapping occurred, we want to go through zero
      if(ffm <= ffp) {
        get_nearby_from_range(low, high, nearby, pos, hh);
      }
      else {
        get_nearby_from_range(sort_los.begin(), high, nearby, pos, hh);
        get_nearby_from_range(low, sort_los.end(), nearby, pos, hh);
      }
}

/*This function takes a particle position and returns a list of the indices of lines near it in index_nr_lines
 * Near is defined as: dx^2+dy^2 < 4h^2 */
std::map<int,double> IndexTable::get_near_lines(const float pos[],const double hh)
{
      std::map<int, double> nearby;
      if(index_table.size() > 0){
        get_nearby(pos[0],index_table, nearby, pos, hh);
      }
      if(index_table_xx.size() > 0){
        get_nearby(pos[1],index_table_xx, nearby, pos, hh);
      }
      return nearby;
}
