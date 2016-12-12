/*
 *  C program for the FES prediction software.
 *  Adapted version of CSR2mod to FES atlas
 *  File      : fes2mod.c
 *  Date      : february 2007 
 *  Added by IK 
 *  This prog reads info about the model grid lon-lat-depth-dimensions
 *  and amplitudes and phases for 8 constituents from a global grid,and interpolates
 *  these to the model boundary grid points and save the result to a file
 *  FESobc_elev.dat which is read by HYCOM to compute the tidal boundary conditions.
 *  It reads the following files:
 *   newpos.uf         --- contains model lon lat positions.
 *   depths???x???.uf  --- model bathymetry.
 *   
 * ----------------------------------------------------------------------
 *
 *  rc      : return code (problem if rc != 0)
 *  lat     : latitude
 *  lon     : longitude
 *
 */       

#include <stdlib.h>
#include <stdio.h>
#include "fes.h"
#include "fes-int.h" // in order to use EPSILON

static const int nconst = 8; 
static const int nx     = 255; 
static const int ny     = 225; 

// Order constituents list for csr:
// Q1 O1 P1 K1 N2 M2 S2 K2
static const char* namewave[] =
{
 "Q1  ", "O1  ", "P1  ","K1  ", "N2  ", "M2  ", "S2  ", "K2  "
};

// corresponding wave code number in FES2004 is respectively 
static const int iwave[] = 
{
 0, 1, 8, 2, 4, 5, 7, 6
};


int main(void)
{ 
  fesData*	shortTide = NULL;
  char* fespath;
  int		rc	= 0;
  FILE          *fout;
  int           i,k,iamp,ipha; 
  float        dlon,dlat;
  double        ampha[2*nconst];

// Get the frequency in the same units than in mod_tides.F90 on hycom: cycles/day
  double freq[ ] =
  {
   0.893244048,0.929535695,0.997262079,1.002737898,1.895981946,1.932273593,2.000000000,2.005475795
  };


 
  printf("Enter lat and lon\n");
  scanf("%f %f",&dlat,&dlon);
  //scanf("%d %d",&i,&k);
  if ( (abs(dlat) < EPSILON) | (abs(dlon) < EPSILON) )
  {
     printf("Try again and Enter lat lon !!! %f %f \n",dlat,dlat);
     //printf("Try again and Enter lat lon !!! %d %d \n",i,k);
     return -1;
  }



 /* KAL - added option to specify FES_PATH */
  if ( getenv("FES2004_PATH") == NULL ) {
     printf("\nError: Set enironment variable FES2004_PATH to the location of FES2004 data \n");
     return -1;
  } else {
     fespath=getenv("FES2004_PATH") ;
  }


  if ( rc != FES_SUCCESS )
      goto onError;


 /*--------- Initialize memory for FES algorithms ------------------*/
  //rc = fesNew(&shortTide, FES_TIDE, FES_IO, "/home/fimm/nersc/intissar/Progs/Tide_FES/data");
  rc = fesNew(&shortTide, FES_TIDE, FES_IO, fespath);
  if ( rc != FES_SUCCESS )
      goto onError;


/*-------------------------- Estimate boundary values -------------------*/
/* NB: in C/C++ arrays start at 0!!!                                   --*/
/* --- Data are generated mlong for k=1, and k=nx-1 mlong the          --*/
/* --- 'EW' boundaries and mlong j=1, j=jj-2 mlong the 'NS' boundaries.--*/
/* --- If this is changed, the calculations on index m below, and in   --*/
/* --- subroutine nest_barotp.F90  in the model code must be modeified.--*/


  rc=diag8Const(shortTide,dlat,dlon,nconst,&ampha);
     if ( rc != FES_SUCCESS )
        goto onError;

/**************************   DUMP DATA IN FESobc_elev.dat  *********************/


 fout=fopen("testlatlon.dat","w");

 fprintf(fout,"'Name Freq (cyc/d)'\n");
 for (i=0;i<nconst;i++)
 {
    k=iwave[i];
    fprintf(fout,"'%s %10.9f'\n",namewave[i],freq[i]);
 }
 fprintf(fout,"    Lat    Lon       Amp    Phas\n");
 for (i=0;i<nconst;i++) 
 {
    iamp=2*i;
    ipha=2*i+1;
    fprintf(fout,"%9.4f %9.4f %9.4f %9.4f\n",dlat,dlon,ampha[iamp],ampha[ipha]);
 }
 printf("Data dumped in file testlatlon.dat\n"); 
 

 goto onTerminate;

onError:
  printf("#ERROR : %s\n", fesErr2String(rc)),
  rc = 1;

onTerminate:
  /* Free memory for FES algorithms */
  fesDelete(shortTide);
  fclose(fout);
  return rc;
}
