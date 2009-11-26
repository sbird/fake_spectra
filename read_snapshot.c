
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
  char input_fname[200];
  int type, files;
  
  type = PARTTYPE; 
  files = FILENUMBER;
  
  if(argv[1] == NULL)
    {
      printf("Please enter the name of the file you want to analyse.\n");
      printf("exiting.......\n");
      exit(99);
    }
  
  strcpy(input_fname, argv[1]);
  
  load_snapshot(input_fname, files, type);
  
  if(type == 0)
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
  float     *temp;
  int        i,k,ntot_withmasses;
  int        n;
  int Nstart=0;
  
  for(i=0; i<files; i++)
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
      
     if(!read_gadget_head(&header1,fd,0))
     {
        fprintf(stderr, "Can't read header!\n");
        exit(1);
     }
      
      /*Load in the number of different particle types*/
      for(k=0, NumPartTot=0, ntot_withmasses=0; k<5; k++)
	{
	  NumPart[k] = header1.npart[k]; /* Total particle number */
	  NumPartTot += header1.npart[k];
	}
      
      if(NumPart[itype] == 0)
	{
	  printf("File %s contains no particles of type  %d! Exiting.\n\n", buf,itype);
	  exit(0);
	}
      
      for(k=0, ntot_withmasses=0; k<6; k++)
	{
	  if(header1.mass[k]==0)
	    ntot_withmasses+= header1.npart[k];
	} 
      
      for(k=0 ;k<5 ;k++)
	printf("Massarr[%d] = %e\n",k,header1.mass[k]);
      printf("Particles with mass = %d\n\n",ntot_withmasses);
      
      
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
      for(n=0; n<itype; n++)
              Nstart+=NumPart[n];
      if(i==0)
	allocate_memory();
  
      if(!(temp=malloc(3*Ntype*sizeof(float))))
      {
                      fprintf(stderr, "Failed to allocate temp memory\n");
                      exit(2);
      }
      /*Particle positions */
      printf("Reading position...checking\n");
      read_gadget_float3(temp,"POS ",Nstart,Ntype,fd,0);
      for(n=0;n<Ntype;n++)
      {
         P[n].Pos[0]=temp[3*n];	   
         P[n].Pos[1]=temp[3*n+1];	   
         P[n].Pos[2]=temp[3*n+2];	   
      }
      
      printf("P[%d].Pos[0] = %f\n", 0, P[0].Pos[0]);
      printf("P[%d].Pos[1] = %f\n", 0, P[0].Pos[1]);
      printf("P[%d].Pos[2] = %f\n\n", 0, P[0].Pos[2]);
      
      /* Peculiar velocites */
      printf("Reading velocity...checking\n");
      read_gadget_float3(temp,"VEL ",Nstart,Ntype,fd,0);
      for(n=0;n<Ntype;n++)
      {
         P[n].Vel[0]=temp[3*n];	   
         P[n].Vel[1]=temp[3*n+1];	   
         P[n].Vel[2]=temp[3*n+2];	   
      }
      printf("P[%d].Vel[0] = %f\n", 0, P[0].Vel[0]);
      printf("P[%d].Vel[1] = %f\n", 0, P[0].Vel[1]);
      printf("P[%d].Vel[2] = %f\n\n", 0, P[0].Vel[2]);
      free(temp);
      if(!(temp=malloc(Ntype*sizeof(float))))
      {
                      fprintf(stderr, "Failed to allocate temp memory\n");
                      exit(2);
      }
      
      /* Particles masses  */
      printf("Reading mass...checking\n");
      if(ntot_withmasses>0)
              read_gadget_float(temp,"MASS",Nstart,Ntype,fd);
      for(n=0;n<Ntype;n++)
      {
   	  if(header1.mass[itype]>0.0)
   	      P[n].Mass = header1.mass[itype];
   	  else				
   	      P[n].Mass=temp[n];
      }
      
      printf("P[%d].Mass = %e\n\n", Ntype, P[1].Mass);
      
      
      if(itype == 0)
	{ 
	  /*The internal energy of all the Sph particles is read in */
	  printf("Reading internal energy per unit mass...checking\n");
              read_gadget_float(temp,"U   ",0,NumPart[0],fd);
	  for(n=0; n<NumPart[0];n++)
	      P[n].U=temp[n];
	  
	  printf("P[%d].U = %f\n\n", Ntype, P[1].U);
	  
          /* The free electron fraction */
          if(header1.flag_cooling)
            {
              printf("Reading electron fraction...checking\n");
              read_gadget_float(temp,"NHP ",0,NumPart[0],fd);
	      for(n=0; n<NumPart[0];n++)
	         P[n].Ne=temp[n];

              printf("P[%d].Ne = %e\n\n", Ntype, P[1].Ne);
              /* The HI fraction, nHI/nH */
              read_gadget_float(temp,"NH  ",0,NumPart[0],fd);
              for(n=0; n<NumPart[0];n++)
	        P[n].NH0=temp[n];
	      printf("P[%d].NH0 = %e\n\n", Ntype, P[1].NH0);
	    }
 
	  /* The smoothing length */	  
         read_gadget_float(temp,"HSML",0,NumPart[0],fd);
	 for(n=0; n<NumPart[0];n++)
	    P[n].h=temp[n];
	 printf("P[%d].h = %f\n\n",Ntype, P[1].h);
	  
	}
      fclose(fd);
      free(temp);
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
  
  printf("Allocating memory...done\n");
  return(0);
}
/*****************************************************************************/
