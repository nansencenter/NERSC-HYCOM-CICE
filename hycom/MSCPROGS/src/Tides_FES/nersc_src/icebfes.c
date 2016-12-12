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

#include "fes.h"

int main(void)
{
  fesData*	shortTide = NULL;
  fesData*	radialTide = NULL;
  char* fespath;
  int		rc	= 0;
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
     
  /* Read the data file first */

  fdata=fopen("fes_ptt3106.txt","r");
  fout=fopen("height_ptt3106.txt","w");
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
      printf("Could not open the file");
      goto onTerminate;
  }

 /* KAL - added option to specify FES_PATH */
  if ( getenv("FES2004_PATH") == NULL ) {
     printf("\nError: Set enironment variable FES2004_PATH to the location of FES2004 data \n");
     return -1;
  } else {
     fespath=getenv("FES2004_PATH") ;
  }

  npos=k;
  printf("nb position %d,%f,%f\n",k,time[k],lat[k]);
  printf("\n\n");
  /* Initialize memory for FES algorithms */
  rc = fesNew(&shortTide, FES_TIDE, FES_IO, fespath);
  if ( rc != FES_SUCCESS )
      goto onError;

  rc = fesNew(&radialTide, FES_RADIAL, FES_IO, fespath);
  if ( rc != FES_SUCCESS )
      goto onError;

 /* printf("%12s %5s %9s %9s %9s %9s %9s %9s %9s\n",
    "JulDay","Hour","Latitude","Longitude",
    "Short_tid","LP_tid","Pure_Tide","Geo_Tide","Rad_Tide");*/
  //for(k = 1; k <npos ; k++)
  for(k = 0; k <npos ; k++)
  {
    /* Compute tide */
    rc = fesCore(shortTide, lat[k], lon[k], time[k], &tide, &lp);
    if ( rc == FES_SUCCESS )
    {
      /* Compute load tide */
      rc = fesCore(radialTide, lat[k], lon[k], time[k], &load, &loadlp);
      if(  rc != FES_NO_DATA && rc == FES_SUCCESS )
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
    else
    {
      goto onError;
    }
  }

  goto onTerminate;

onError:
  printf("#ERROR : %s\n", fesErr2String(rc)),
  rc = 1;

onTerminate:
  /* Free memory for FES algorithms */
  fesDelete(shortTide);
  fesDelete(radialTide);
  fclose(fdata);
  fclose(fout);

  return rc;
}
