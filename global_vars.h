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
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];
} header1;


int    NumPart[6], Ntype, NumPartTot;

struct particle_data 
{
  float  Pos[3];
  float  Vel[3];
  float  Mass;
  int    Type;
  double Mass_d;
  float  Rho, U, NH0, NHep, Ne, h;
} *P;

int *Id;

double  atime, redshift, omega0, omegaL, box100, h100;
