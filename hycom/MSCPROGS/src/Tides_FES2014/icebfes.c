/*
 *  C program for the FES prediction software.
 * 
 *  File      : testfes.c
 *  Developer : CLS - CNRS/LEGOS
 *  Version   : 1.3
 *  Date      : 13 January 2005
 *  
 * ----------------------------------------------------------------------
 *
 *  rc      : return code (problem if rc != 0)
 *  lat     : latitude
 *  lon     : longitude
 *  time    : time in CNES Julian days
 *  hour    : hour
 *  tide    : short tides (semi_diurnal and diurnal tides)
 *  lp      : long period tides
 *  load    : loading effects for short tide
 *  loadlp  : loading effects for long period tides
              (is always equal to zero)
 * 
 *  tide+lp             = pure tide (as seen by a tide gauge)
 *  tide+lp+load        = geocentric tide ((as seen by a satellite)
 *  CNES Julian day     = NASA Julian day + 2922
 *  CNES Julian day 0 is at midnight between the 31 December 
 *                          and 01 January 1950 AD Gregorian
 *
 */       

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "fes.h"

int main(void)
{
  FES	shortTide = NULL;
  FES	radialTide = NULL;
  char* fespath;
  double	tide;
  double	lp;
  double	load;
  double	loadlp;

  float  	lon[500];
  float	        lat[500];
  float  	time[500];	/* 1983-01-01 00:00:00.0 */

  int		npos	= 0;
  int           k       = 0;
  FILE*        fdata = NULL;
  FILE*        fout = NULL;
  static const char* filein="fes_ptt3106.txt";
  static const char* fileout="height_ptt3106.txt";
  int errtide, errload;
     
  /* Read the data file first */
  fdata=fopen(filein,"r");
  fout=fopen(fileout,"w");
  if (fdata != NULL)
  {  
     while (fscanf(fdata,"%f %f %f \n",&time[k],&lat[k],&lon[k]) != EOF)
     //while (fscanf(fdata,"%f %f %f \n",&t,&la,&lo) != EOF)
     {
        printf("line number: %d %f %f %f\n",k,time[k],lat[k],lon[k]);
        k=k+1;
     }
 
  }
  else
  {
      printf("Could not open the file %s\n",filein);
      exit(1);
  }
  fclose(fdata);

 /* KAL - added option to specify FES_PATH */
  if ( getenv("FES2014_PATH") == NULL ) {
     printf("\nError: Set enironment variable FES2014_PATH to the location of FES2014 data \n");
     return -1;
  } else {
     fespath=getenv("FES2014_PATH") ;
  }
  printf("FES2014_PATH set to %s\n",fespath);
  setenv("FES_DATA",fespath,1);
  char * tmp = "/ocean_tide.ini";
  char * tmp2= "/load_tide.ini";
  int newSize = strlen(fespath)  + strlen(tmp) + 1; 
  char * inifile     = (char *)malloc(newSize);
  newSize = strlen(fespath)  + strlen(tmp2) + 1; 
  char * inifileload = (char *)malloc(newSize);


  strcpy(inifile,fespath);
  strcat(inifile,tmp); // or strncat 
  printf("FES ini file is set to %s\n",inifile);
  strcpy(inifileload,fespath);
  strcat(inifileload,tmp2); // or strncat 
  printf("FES radial ini file is set to %s\n",inifileload);

  npos=k;
  /* Initialize memory for FES algorithms */
  if (fes_new(&shortTide, FES_TIDE, FES_IO, inifile) )
  {
     printf("fes_new error : %s\n", fes_error(shortTide));
     return 1;
  }

  if ( fes_new(&radialTide, FES_RADIAL, FES_IO, inifileload))
  {
     printf("fes_new error : %s\n", fes_error(radialTide));
     return 1;
  }

 /* printf("%12s %5s %9s %9s %9s %9s %9s %9s %9s\n",
    "JulDay","Hour","Latitude","Longitude",
    "Short_tid","LP_tid","Pure_Tide","Geo_Tide","Rad_Tide");*/
  //for(k = 1; k <npos ; k++)
  for(k = 0; k <npos ; k++)
  {

    errtide=0;
    errload=0;


    if ((errtide=fes_core(shortTide, lat[k], lon[k], time[k], &tide, &lp)))
    {
       printf("fes_new error tide  : %s\n", fes_error(shortTide));
    }

    if ((errload=fes_core(radialTide, lat[k], lon[k], time[k], &load, &loadlp)))
    {
       printf("fes_new error load  : %s\n", fes_error(radialTide));
    }

    if (errtide==0 && errload==0) 
    {
       fprintf(fout,"%12.5f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n",
         time[k],
         lat[k],
         lon[k], 
         tide,
         lp,
         tide+lp,
         tide+lp+load,
         load);
   }
   }
   fclose(fout);
   //fprintf("Results in %s\n",fileout);


   /* Free memory for FES algorithms */
   fes_delete(shortTide);
   fes_delete(radialTide);
   return 0;
}
