#ifndef INDEX_TABLE_H
#define INDEX_TABLE_H

#include <map>
#include "types.h"

class IndexTable
{
public:
  IndexTable(const los *los_table_i, const int NumLos_i, const double box);
  //Get a list of lines nearby a particle with coordinates xx, yy, zz and smoothing length hh.
  std::map<int,double> get_near_lines(const float pos[], const double hh);
private:
  //Get a list of lines nearby a particle from a particular index table
  void get_nearby(float first, std::map<double, int>& sort_los, std::map<int, double>& nearby, const float pos[], const double hh);
  //Get a list of lines nearby a particle from an iterator range
  void get_nearby_from_range(std::map<double, int>::iterator low, std::map<double, int>::iterator high, std::map<int, double>& nearby, const float pos[], const double hh);

  // The key is the position of the primary axis, which is xx for index_table and yy for index_table_xx.
  // The value is the index of this entry in los_table.
  //index_table stores lines where axis = 2 or 3.
  std::map<double, int> index_table;
  //index_table_xx stores lines where axis = 1.
  std::map<double, int> index_table_xx;
  //Pointers to the original los table
  const los *los_table;
  const int NumLos;
  const double boxsize;
};

void Compute_Absorption(double * tau_H1, double * rho, double * veloc, double * temp, const int nbins, const double Hz, const double h100, const double box100, const double atime, const double lambda_lya, const double gamma_lya, const double fosc_lya, const double mass);

void SPH_Interpolation(double * rhoker_H, interp * species, const int nspecies, const int nbins, const int Particles, const int NumLos,const double boxsize, const los *los_table, IndexTable& sort_los_table, const pdata *P);
void Rescale_Units(double * rho, double * veloc, double * temp, const int nbins, const double h100, const double atime);
void Convert_Density(double * rhoker_H, double * rho, const double h100, const double atime, const double omegab);

#endif
