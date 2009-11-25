
#include "headers.h"
#include "global_vars.h"
#include "parameters.h"

int allocate_memory();
int load_snapshot();
void SPH_interpolation();

/* Here we load a snapshot file. It can be distributed
 * onto several files (for files>1)
 */
/*****************************************************************************/
int main(int argc, char **argv)
{
  char path[200], input_fname[200], basename[200];
  int  type, snapshot_number, files;
  int i;
  
  itype = PARTTYPE; 
  files = FILENUMBER;
  
  if(argv[1] == NULL)
    {
      printf("Please enter the name of the file you want to analyse.\n");
      printf("exiting.......\n");
      exit(99);
    }
  
  strcpy(input_fname, argv[1]);
  
  load_snapshot(input_fname, files, itype);
  
  if(itype == 0)
    SPH_interpolation();
  
  return(0);
}
/*****************************************************************************/
/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 */
int load_snapshot(char *fname, int files, int itype)
{
  FILE *fd;
  char   buf[200];
  int    i,j,k,dummy,ntot_withmasses;
  int    t,n,off,pc,pc_new,pc_sph;
  int ii=0;
  float dummy2;
  
  
  
#define SKIP fread(&dummy, sizeof(dummy), 1, fd);
  
  for(i=0, pc=1; i<files; i++, pc=pc_new)
    {
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
      
      
      /* Read the header data */
      fread(&dummy, sizeof(dummy), 1, fd);
      fread(&header1, sizeof(header1), 1, fd); 
      fread(&dummy, sizeof(dummy), 1, fd);
      
      
      /*Load in the number of different particle types*/
      for(k=0, NumPartTot=0, ntot_withmasses=0; k<5; k++)
	{
	  NumPart[k] = header1.npart[k]; /* Total particle number */
	  NumPartTot += header1.npart[k];
	}
      
      if(NumPart[itype] == 0)
	{
	  printf("File %s contains no particles! Exiting.\n\n", buf);
	  exit(0);
	}
      
      for(k=0, ntot_withmasses=0; k<5; k++)
	{
	  if(header1.mass[k]==0)
	    ntot_withmasses+= header1.npart[k];
	} 
      
      for(k=0 ;k<5 ;k++)
	printf("Massarr[%d] = %e\n",k,header1.mass[k]);
      
      
      atime= header1.time;
      redshift= header1.redshift;
      box100 = header1.BoxSize;
      omega0 = header1.Omega0;
      omegaL = header1.OmegaLambda;
      h100 = header1.HubbleParam;
      
    
      printf("NGas = %d\n",NumPart[0]);
      printf("NDM = %d\n", NumPart[1]);
      printf("NStars = %d\n", NumPart[4]);
      printf("Total particle number = %d\n\n",NumPartTot);
      printf("Redshift = %f\n",redshift);
      printf("Expansion factor = %f\n",atime);
      printf("Lambda = %f\n",omegaL);
      printf("Matter = %f\n",omega0);
      printf("Hubble = %f\n",h100);
      printf("Box = %f\n\n",box100);
    
      
      /* Now read the particle data */ 
      Ntype = NumPart[itype];
      
      if(i==0)
	allocate_memory();
      
      /*Particle positions */
      printf("Reading position...checking\n");
      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<NumPart[k];n++)
	    {
	      if(k == itype)						
		{
		  fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);	   
		  pc_new++;
		}
	      else
		for (ii=0;ii<3;ii++)
		  fread(&dummy2, sizeof(dummy2), 1, fd);	      
	    }
	  
	}
      SKIP;
      
      printf("P[%d].Pos[0] = %f\n", 1, P[1].Pos[0]);
      printf("P[%d].Pos[1] = %f\n", 1, P[1].Pos[1]);
      printf("P[%d].Pos[2] = %f\n\n", 1, P[1].Pos[2]);
      
      /* Peculiar velocites */
      printf("Reading velocity...checking\n");
      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<NumPart[k];n++)
	    {
	      if(k == itype)						
		{
		  fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);
		  pc_new++;
		}
	      else
		for (ii=0;ii<3;ii++)
		  fread(&dummy2, sizeof(dummy2), 1, fd);	      
	    }
	}
      SKIP;

      printf("P[%d].Vel[0] = %f\n", 1, P[1].Vel[0]);
      printf("P[%d].Vel[1] = %f\n", 1, P[1].Vel[1]);
      printf("P[%d].Vel[2] = %f\n\n", 1, P[1].Vel[2]);
      
      
      /* Particle IDs (unused)  */
      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	for(n=0;n<NumPart[k];n++)
	  fread(&dummy, sizeof(dummy), 1, fd);
       SKIP;
      
      
      /* Particles masses  */
      printf("Reading mass...checking\n");
      if(ntot_withmasses>0)
       	SKIP;
      for(k=0, pc_new=pc; k<6; k++)
	{
	  for(n=0;n<NumPart[k];n++)
	    {
	      if(k == itype)				
		{
		  if(header1.mass[k]>0.0)
		    {
		      P[pc_new].Mass = header1.mass[k];
		      pc_new++;			
		    }
		  else				
		    {
		      fread(&P[pc_new].Mass, sizeof(float), 1, fd);
		      pc_new++;
		    }
		}
	      else
		{
		  if(header1.mass[k]>0.0)
		    dummy2 = header1.mass[k];
		  else
		    fread(&dummy2, sizeof(dummy2), 1, fd);
		}
	      
	    }
	}
      if(ntot_withmasses>0)
	SKIP;
      
      printf("P[%d].Mass = %e\n\n", Ntype, P[1].Mass);
      
      
      if(itype == 0)
	{ 
	   SKIP;
	  /*The internal energy of all the Sph particles is read in */
	  printf("Reading internal energy per unit mass...checking\n");
	  for(n=0, pc_sph=pc; n<NumPart[0];n++)
	    {
	      fread(&P[pc_sph].U, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP;
	  
	  printf("P[%d].U = %f\n\n", Ntype, P[1].U);
	  
	  
	  
	 SKIP;

	  /*The comoving density of all the SPH particles is read in */
          printf("Reading comoving SPH density...checking\n");
          for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
              fread(&P[pc_sph].Rho, sizeof(float), 1, fd);

              pc_sph++;
            }

	  printf("P[%d].Rho = %e\n", Ntype, P[1].Rho);
          printf("\n");
	  SKIP;
	  

          /* The free electron fraction */
          if(header1.flag_cooling)
            {
	      SKIP;
              printf("Reading electron fraction...checking\n");
              for(n=0, pc_sph=pc; n<NumPart[0];n++)
                {
                  fread(&P[pc_sph].Ne, sizeof(float), 1, fd);
                  pc_sph++;
                }
              SKIP;
            }

          printf("P[%d].Ne = %e\n\n", Ntype, P[1].Ne);



      	  /* The HI fraction, nHI/nH */
	  if(header1.flag_cooling)
	    {
	      SKIP;
	      
	      printf("Reading HI fraction...checking\n");
	      for(n=0, pc_sph=pc; n<NumPart[0];n++)
		{
		  fread(&P[pc_sph].NH0, sizeof(float), 1, fd);
		  pc_sph++;
		}
	        SKIP;
	    }
	  
	  printf("P[%d].NH0 = %e\n\n", Ntype, P[1].NH0);
	  
 
	  /* The smoothing length */	  
	   SKIP;
	  printf("Reading smoothing length...checking\n");
	  for(n=0, pc_sph=pc; n<NumPart[0];n++)
	    {
	      fread(&P[pc_sph].h, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP;
	  
	  printf("P[%d].h = %f\n\n",Ntype, P[1].h);
	  
	  
	  
	}
      fclose(fd);
    }
      

  printf("Baryon fraction in box = %f\n\n",OMEGAB);
  
 
  return(0);
}
/*****************************************************************************/
/* this routine allocates the memory for the 
 * particle data.
 */
  int allocate_memory(void)
  { 
  if(!(P=malloc(Ntype*sizeof(struct particle_data))))
    {
      fprintf(stderr,"failed to allocate memory.\n\n");
      exit(0);
    }
  
  P--;   /* start with offset 1 */

  
  if(!(Id=malloc(Ntype*sizeof(int))))
    {
      fprintf(stderr,"failed to allocate memory.\n\n");
      exit(0);
    }
  
  Id--;   /* start with offset 1 */

  printf("Allocating memory...done\n");
  return(0);
}
/*****************************************************************************/
