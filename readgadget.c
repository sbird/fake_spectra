#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>
#include <stdint.h>
#include <sys/types.h>
#include <unistd.h>
#include "global_vars.h"
#define int4bytes int
gadget_header header;

/*--------- comment/uncomment to remove/enable DEBUG outputs ------------------*/
/*
#define MY_DEBUG
*/

/*-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/
/*---------------------------- Low Level Routines -----------------------------*/
/*-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/

int4bytes blksize,swap=0;
int nparticles=0;
#define SKIP  {my_fread(&blksize,sizeof(int),1,fd); swap_Nbyte((char*)&blksize,1,4);}

/*---------------------- Basic routine to read data from a file ---------------*/
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nread;

  if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
    {
      int err;
      if(feof(stream))
              return 0;
      fprintf(stderr, "fread error: %zd = fread(%d %zd %zd file)!\n",nread,fileno(stream),size,nmemb);
      err=ferror(stream);
      fprintf(stderr, "ferror gives: %d : %s\n",err,strerror(err));
      exit(3);
    }
  return nread;
}

/*-----------------------------------------------------------------------------*/
/*---------------------- Routine to swap ENDIAN -------------------------------*/
/*-------- char *data:    Pointer to the data ---------------------------------*/
/*-------- int n:         Number of elements to swap --------------------------*/
/*-------- int m:         Size of single element to swap ----------------------*/
/*--------                int,float = 4 ---------------------------------------*/
/*--------                double    = 8 ---------------------------------------*/
/*-----------------------------------------------------------------------------*/
void swap_Nbyte(char *data,int n,int m)
{
  int i,j;
  char old_data[16];

  if(swap>0)
    {
      for(j=0;j<n;j++)
	{
          memcpy(&old_data[0],&data[j*m],m);
          for(i=0;i<m;i++)
            {
              data[j*m+i]=old_data[m-i-1];
	    }
	}
    }
}

/*-----------------------------------------------------------------------------*/
/*---------------------- Routine find a block in a snapfile -------------------*/
/*-------- FILE *fd:      File handle -----------------------------------------*/
/*-------- char *label:   4 byte identifyer for block -------------------------*/
/*-------- returns length of block found, -------------------------------------*/
/*-------- the file fd points to starting point of block ----------------------*/
/*-----------------------------------------------------------------------------*/
int_blk find_block(FILE *fd,char *label)
{
  int_blk blocksize=0;
  char blocklabel[5]={"    "};

  rewind(fd);

  while(blocksize == 0)
    {
       SKIP;
       if(blksize == 134217728)
	 {
#ifdef MY_DEBUG
           printf("Enable ENDIAN swapping !\n");
#endif
           swap=1-swap;
           swap_Nbyte((char*)&blksize,1,4);
	 }
       if(blksize != 8)
         {
	   fprintf(stderr,"incorrect format (blksize=%d)!\n",blksize);
           exit(1);
         }
       else
         {
           my_fread(blocklabel, 4*sizeof(char), 1, fd);
           if(feof(fd)) 
           {
             fprintf(stderr, "Block (float) <%s> not found!\n",label);
             exit(5);
           }
           my_fread(&blocksize, sizeof(int4bytes), 1, fd);
           if(blocksize > pow(2,32))
           {
              fprintf(stderr, "You are apparently trying to read a very large file. This is not (yet) supported.\n");
              exit(1);
           }
           swap_Nbyte((char*)&blocksize,1,4);
#ifdef MY_DEBUG
         fprintf(stderr, "Found Block <%s> with %ld bytes\n",blocklabel,blocksize);
#endif
           SKIP;
	   if(strcmp(label,blocklabel)!=0)
	     { 
                fseek(fd,blocksize,1);
                blocksize=0;
	     }
         }
    }
  return(blocksize-8);
}


/*-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/
/*---------------------------- High Level Routines ----------------------------*/
/*-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------*/
/*---------------------- Routine to read the header information ---------------*/
/*-------- int *npart:    List of Particle numbers for spezies [0..5] ---------*/
/*-------- int *massarr:  List of masses for spezies [0..5] -------------------*/
/*-------- int *time:     Time of snapshot ------------------------------------*/
/*-------- int *redshift:  Redshift of snapshot --------------------------------*/
/*-------- FILE *fd:      File handle -----------------------------------------*/
/*-------- returns number of read bytes ---------------------------------------*/
/*-----------------------------------------------------------------------------*/
int read_gadget_head(gadget_header *out_header, FILE *fd, int old)
{
  int_blk blocksize;
  if(!old)
    blocksize = find_block(fd,"HEAD");
  else{
    blocksize=sizeof(header);
    blksize=8;
  }
  if(blocksize <= 0)
    {
      fprintf(stderr, "Block <%s> not found!\n","HEAD");
      exit(5);
    }
  else
    {
				SKIP;
		 my_fread(&header, sizeof(header), 1, fd);
	#ifdef MY_DEBUG
		 fprintf(stderr, "npart=%d %d %d %d %d %d\n", header.npart[0],header.npart[1],header.npart[2],header.npart[3], header.npart[4], header.npart[5]);
		 fprintf(stderr, "header.mass=%f %f %f %f %f %f\n", header.mass[0],header.mass[1],header.mass[2],header.mass[3], header.mass[4], header.mass[5]);
		 fprintf(stderr, "time=%f\n",header.time);
		 fprintf(stderr, "redshift=%g\n",header.redshift);
		 fprintf(stderr, "flag_sfr=%d\n",header.flag_sfr);
		 fprintf(stderr, "flag_feedback=%d\n",header.flag_feedback);
		 fprintf(stderr, "header.npartTotal=%d %d %d %d %d %d\n", header.npartTotal[0],header.npartTotal[1],header.npartTotal[2],header.npartTotal[3], header.npartTotal[4], header.npartTotal[5]);
		 fprintf(stderr, "flag_cooling=%d\n",header.flag_cooling);
		 fprintf(stderr, "numfiles=%d\n",header.num_files);
		 fprintf(stderr, "boxsize=%e\n",header.BoxSize);
	#endif
       SKIP;
       *out_header=header;
    }
  return(blocksize);
}

/*-----------------------------------------------------------------------------*/
/*---------------------- Routine to read a 1D float array ---------------------*/
/*-------- float *data:     Pointer where the data are stored to ----------------*/
/*-------- char *label:   Identifyer for the datafield to read ----------------*/
/*-------- FILE *fd:      File handle -----------------------------------------*/
/*-------- returns length of dataarray ----------------------------------------*/
/*-----------------------------------------------------------------------------*/
int_blk read_gadget_float(float *data,char *label,int offset, int number, FILE *fd, int old)
{
  int_blk blocksize;
  if(!old)
    blocksize = find_block(fd,label);
  else{
	blocksize = header.npart[1]*sizeof(float);
  }
  if(blocksize <= 0)
    {
      fprintf(stderr, "Block (float) <%s> not found!\n",label);
      exit(5);
    }
  else
    {
       blocksize=(blocksize < number*sizeof(float) ? blocksize : number*sizeof(float));
#ifdef MY_DEBUG
       printf("Reading %ld bytes of data from <%s>...\n",blocksize,label);
#endif
       SKIP;
       if(offset>0)
          fseek(fd,offset*sizeof(float),SEEK_CUR);
       my_fread(data,blocksize, 1, fd);
       swap_Nbyte((char*)data,blocksize/sizeof(float),4);
       SKIP;
    }
  return(blocksize/sizeof(float));
}

/*-----------------------------------------------------------------------------*/
/*---------------------- Routine to read a 3D float array ---------------------*/
/*-------- float *data:     Pointer where the data are stored to ----------------*/
/*-------- char *label:   Identifyer for the datafield to read ----------------*/
/*-------- int offset:   Number of elements into the block to start. 
 *-------- Note an "element" is 3 floats. ----------------*/
/*-------- int number:  Total values to read ----------------*/
/*-------- FILE *fd:      File handle -----------------------------------------*/
/*-------- returns length of dataarray ----------------------------------------*/
/*-----------------------------------------------------------------------------*/
int_blk read_gadget_float3(float *data,char *label,int offset, int number, FILE *fd, int old)
{
  int_blk blocksize;

  if(!old)
    blocksize = find_block(fd,label);
  else{
	blocksize=3*header.npart[1]*sizeof(float);
  }
  if(blocksize <= 0)
    {
      fprintf(stderr, "Block (float3) <%s> not found!\n",label);
      exit(5);
    }
  else
    {
       blocksize=(blocksize < 3*number*sizeof(float) ? blocksize : 3*number*sizeof(float));
#ifdef MY_DEBUG
       fprintf(stderr,"Reading %ld bytes of data from <%s>...\n",blocksize,label);
#endif
       SKIP;
       if(offset>0)
          fseek(fd,3*offset*sizeof(float),SEEK_CUR);
       my_fread(data,blocksize, 1, fd);
       swap_Nbyte((char*)data,blocksize/sizeof(float),4);
#ifdef MY_DEBUG
		 fprintf(stderr, "first particles at: %e %e %e\n",data[0],data[1],data[2]);
		 fprintf(stderr, "last particles at: %e %e %e\n",data[number-3],data[number-2],data[number-1]);
#endif
       SKIP;
    }
  return(blocksize/sizeof(float)/3);
}

/*-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/
/*---------------------------- Small Example HowToUse -------------------------*/
/*-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/

#ifdef MAKE_READGADGET_TEST
int main(int argc, char **argv)
{
  FILE *fd = 0;
  char filename[265];
  int i,n,ntot;
  int npart[6];
  double masstab[6],redshift,time;
  float *rho;
  struct vector *pos,*b;

  if(argc >= 2)
    {
      strcpy(filename,argv[1]);
      if(!(fd = fopen(filename,"r")))
        {
           printf("Cant open file <%s> !\n",filename);
           exit(2);
        }  
      else
        {
	  /*----------- READ HEADER TO GET GLOBAL PROPERTIES -------------*/
 	   n = read_gadget_head(npart,masstab,&time,&redshift,fd);

           ntot=0;
           for(i=0;i<6;i++)
	     {
	       printf("PartSpezies %d, anz=%d, masstab=%f\n",i,npart[i],masstab[i]);
               ntot += npart[i];
	     }
           printf("Time of snapshot=%f, z=%f, ntot=%d\n",time,redshift,ntot);

	   /*---------- ALLOCATE MEMORY ---------------------------------*/
           rho=malloc(npart[0]*sizeof(float));
           b=malloc(3*npart[0]*sizeof(float));
           pos=malloc(3*ntot*sizeof(float));

	   /*---------- READ DATA BLOCKS --------------------------------*/
	   n = read_gadget_float(rho,"RHO ",fd);
	   n = read_gadget_float3((float*)pos,"POS ",fd);
	   n = read_gadget_float3((float*)b,"BFLD",fd);
	   /*
           for(i=0;i<npart[0];i++)
	     {
	        printf("%d: (%f %f %f) =  %f\n",i,pos[i].x,pos[i].y,pos[i].z,rho[i]);
	     }
	   */
	   /*---------- CLOSE FILE AND FREE DATA ------------------------*/
	   fclose(fd);

           free(rho);
           free(b);
           free(pos);
	}
    }
  else
    {
      printf("Please give a filename ...\n");
      exit(4);
    }
  exit(0);
} 
#endif










