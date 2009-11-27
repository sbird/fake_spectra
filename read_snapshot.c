
#include "headers.h"
#include "global_vars.h"
#include "parameters.h"

int load_snapshot(char *fname, int files);
void SPH_interpolation(int NumLos, int Ntype);

/* Here we load a snapshot file. It can be distributed
 * onto several files (for files>1) */
/*****************************************************************************/
int main(int argc, char **argv)
{
  int Npart, NumLos, files;
  if(argc<4)
  {
    printf("Usage: ./extract NUMLOS NUMFILES base_filename\n");
    exit(99);
  }
  NumLos=atoi(argv[1]);
  files=atoi(argv[2]);
  if(NumLos <=0 || files <=0)
  {
          printf("Need NUMLOS >0\n");
          exit(99);
  }
  Npart=load_snapshot(argv[3], files);
  if(!PARTTYPE)
    SPH_interpolation(NumLos,Npart);
  return 0;
}
/*****************************************************************************/
/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files. */
int load_snapshot(char *fname, int files)
{
  FILE *fd;
  gadget_header *headers;
  float *temp;
  int i,k,n,NumRead,ntot_withmasses=0;
  int NumPart[6]={0,0,0,0,0,0};
  /*Loop over the headers first, to get totals, then loop over data*/
  headers=malloc(files*sizeof(struct io_header_1));
  if(!headers)
  {
    fprintf(stderr, "Error allocating memory for headers.\n");
    exit(1);
  }
  /*First read all the headers, allocate some memory and work out the totals.*/
  for(i=0; i<files; i++)
  {
    char buf[200];
    if(files>1)
      sprintf(buf,"%s.%d",fname,i);
    else
      sprintf(buf,"%s",fname);
    
    if(!(fd=fopen(buf,"r")))
     {
   		fprintf(stderr,"Error opening file %s for reading!\n", buf);
   		exit(1);
     }
     if(!read_gadget_head(headers+i, fd, 0))
     {
   		fprintf(stderr,"Error reading file header!\n");
   		exit(1);
     }
     //By now we should have the header data.
     fclose(fd);
  }
  atime= headers[0].time;
  redshift= headers[0].redshift;
  box100 = headers[0].BoxSize;
  omega0 = headers[0].Omega0;
  omegaL = headers[0].OmegaLambda;
  h100 = headers[0].HubbleParam;
  /*Check that the headers all match*/
  for(i=0; i<files; i++)
  {
    /*Load in the number of different particle types
     * NumPart[k] holds the total number of particles of that type*/
    for(k=0; k<6; k++)
        NumPart[k] += headers[i].npart[k]; /* Total particle number */
    if(!NumPart[PARTTYPE])
      {
        printf("No particles of type  %d! Exiting.\n\n",PARTTYPE);
        exit(1);
      }
    for(k=0; k<6; k++)
        if(headers[i].mass[k]==0 && headers[i].npart[k])
          ntot_withmasses+= headers[i].npart[k];
    if(atime != headers[i].time ||
       redshift!= headers[i].redshift ||
       box100 != headers[i].BoxSize ||
       omega0 != headers[i].Omega0 ||
       omegaL != headers[i].OmegaLambda ||
       h100 != headers[i].HubbleParam  || 
       headers[i].mass[1] !=headers[0].mass[1] ) 
       {
            fprintf(stderr, "Error. File headers are inconsistent.\n");
            exit(4);
       }
  }
  printf("NumPart=[%d,%d,%d,%d,%d,%d], ",NumPart[0],NumPart[1],NumPart[2],NumPart[3],NumPart[4],NumPart[5]);
  printf("Masses=[%g %g %g %g %g %g], ",headers[0].mass[0],headers[0].mass[1],headers[0].mass[2],headers[0].mass[3],headers[0].mass[4],headers[0].mass[5]);
  printf("Particles with mass = %d\n\n",ntot_withmasses);
  printf("Redshift=%g, Ω_M=%g Ω_L=%g\n",redshift,omega0,omegaL);
  printf("Expansion factor = %f\n",atime);
  printf("Hubble = %g Box=%g \n",h100,box100);

  if(!(P=malloc(NumPart[PARTTYPE]*sizeof(struct particle_data))))
  {
    fprintf(stderr,"failed to allocate memory.\n\n");
    exit(1);
  }
  for(i=0, NumRead=0; i<files; i++)
  {
    char buf[200];
    int Nstart=0;
    int Ntype=0;
    if(files>1)
      sprintf(buf,"%s.%d",fname,i);
    else
      sprintf(buf,"%s",fname);
    
    printf("Reading file : %s\n\n", buf);
    if(!(fd=fopen(buf,"r")))
    {
      printf("can't open file `%s`\n\n",buf);
      exit(0);
    }
    
    /* Now read the particle data */ 
    Ntype = headers[i].npart[PARTTYPE];
    for(n=0; n<PARTTYPE; n++)
            Nstart+=headers[i].npart[n];
    if(!(temp=malloc(3*Ntype*sizeof(float))))
    {
                    fprintf(stderr, "Failed to allocate temp memory\n");
                    exit(2);
    }
    /*Particle positions */
    printf("Reading position and velocity...\n");
    read_gadget_float3(temp,"POS ",Nstart,Ntype,fd,0);
    for(n=0;n<Ntype;n++)
    {
       P[NumRead+n].Pos[0]=temp[3*n];	   
       P[NumRead+n].Pos[1]=temp[3*n+1];	   
       P[NumRead+n].Pos[2]=temp[3*n+2];	   
    }
    
    /* Peculiar velocites */
    read_gadget_float3(temp,"VEL ",Nstart,Ntype,fd,0);
    for(n=0;n<Ntype;n++)
    {
       P[NumRead+n].Vel[0]=temp[3*n];	   
       P[NumRead+n].Vel[1]=temp[3*n+1];	   
       P[NumRead+n].Vel[2]=temp[3*n+2];	   
    }
    free(temp);
    /*Done with triplets, read single valued things*/
    if(!(temp=malloc(Ntype*sizeof(float))))
    {
                    fprintf(stderr, "Failed to allocate temp memory\n");
                    exit(2);
    }
    
    /* Particles masses  */
    printf("Reading mass and temp...\n");
    if(headers[i].mass[PARTTYPE] == 0)
            read_gadget_float(temp,"MASS",Nstart,Ntype,fd);
    for(n=0;n<Ntype;n++)
    {
 	  if(headers[i].mass[PARTTYPE])
 	      P[NumRead+n].Mass = headers[i].mass[PARTTYPE];
 	  else				
 	      P[NumRead+n].Mass=temp[n];
    }
    if(headers[i].mass[0])
          omegab = headers[0].mass[PARTTYPE]/(headers[0].mass[PARTTYPE]+headers[0].mass[1])*omega0;
    else
          omegab = P[0].Mass/(P[0].Mass+headers[0].mass[1])*omega0;
    
    if(PARTTYPE == 0)
      { 
        /*The internal energy of all the Sph particles is read in */
        read_gadget_float(temp,"U   ",0,NumPart[0],fd);
        for(n=0; n<NumPart[0];n++)
            P[NumRead+n].U=temp[n];
        
        printf("P[%d].U = %f\n\n", Ntype, P[1].U);
        
        /* The free electron fraction */
        if(headers[i].flag_cooling)
          {
            printf("Reading electron fractions...\n");
            read_gadget_float(temp,"NHP ",0,NumPart[0],fd);
            for(n=0; n<NumPart[0];n++)
               P[NumRead+n].Ne=temp[n];

            /* The HI fraction, nHI/nH */
            read_gadget_float(temp,"NH  ",0,NumPart[0],fd);
            for(n=0; n<NumPart[0];n++)
              P[NumRead+n].NH0=temp[n];
          }

        /* The smoothing length */	  
       read_gadget_float(temp,"HSML",0,NumPart[0],fd);
       for(n=0; n<NumPart[0];n++)
          P[NumRead+n].h=temp[n];
        
      }
    NumRead+=Ntype;
    fclose(fd);
    free(temp);
  }
  printf("P[%d].Pos = [%g %g %g]\n", 0, P[0].Pos[0], P[0].Pos[1],P[0].Pos[2]);
  printf("P[%d].Vel = [%g %g %g]\n", 0, P[0].Vel[0], P[0].Vel[1],P[0].Vel[2]);
  printf("P[%d].Mass = %e Ω_B=%g\n\n", NumRead, P[1].Mass,omegab);
  printf("P[%d].Ne = %e\n", NumRead, P[1].Ne);
  printf("P[%d].NH0 = %e\n", NumRead, P[1].NH0);
  printf("P[%d].h = %f\n",NumRead, P[1].h);
  return NumRead;
}
