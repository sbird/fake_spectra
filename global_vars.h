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
} header1;


int    NumPart[6], Ntype, NumPartTot;

struct particle_data 
{
  float  Pos[3];
  float  Vel[3];
  float  Mass;
  int    Type;
  double Mass_d;
/*   float  Rho, */
        float U, NH0, NHep, Ne, h;
} *P;


double  atime, redshift, omega0, omegaL, box100, h100;

void swap_Nbyte(char *data,int n,int m);
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream);
int64_t find_block(FILE *fd,char *label);
int64_t read_gadget_float(float *data,char *label,FILE *fd);
/* The final argument, if one, means it will attempt to read an old format file*/
int64_t read_gadget_float3(float *data,char *label,int offset, int read, FILE *fd, int old);
int read_gadget_head(struct io_header_1 *out_header, FILE *fd, int old);
