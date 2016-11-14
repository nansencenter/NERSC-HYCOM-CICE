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
 *  time    : time in CNES Julian days
 *  CNES Julian day     = NASA Julian day + 2922
 *  CNES Julian day 0 is at midnight between the 31 December 
 *                          and 01 January 1950 AD Gregorian
 *
 */       

#include <stdlib.h>
#include <stdio.h>
#include "fes.h"
#include <string.h>


static const int nconst = 8; 
//static const int nx     = 255; 
//static const int ny     = 225; 

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

/*
// ///////////////////////////////////////////////////////////////////////////
// function swaps a "size" byte variable, i.e. it converts their
// representation between little and big endian (in either direction).
//
// Parameters:
//   p:		pointer to the variable to swap
//   n:		size of the variable
*/
static void swapByte(void *p, int n)
{
    register char* begin        = (char*)p;
    register char* end          = ((char*)p) +(n - 1);
    register char  temp;

    while ( begin < end )
    {
	temp = *begin;
	*(begin++) = *end;
	*(end--)   = temp;
    }
}


int max(int a,int b)
{
   if (a>b) return a;
   else return b;
}

int min(int a,int b)
{
   if (a<b) return a;
   else return b;
}

/*
//////////////////////////////////////////////////////////////////////////////

*/
int main(void)
{ 
  fesData*	shortTide = NULL;
  int		rc	= 0;
  FILE          *fout;
  int           nx,ny;
  int           i,j,k,m,itec,i1,i2; 
  double        dlon,dlat;
  FILE*         flatlon=NULL;
  FILE*         fdepths=NULL;
  char          fdepth_name[16];
  char          var[]="VARIABLES=\"m\""; 
  char          line[100];

/* Constituents frequency defined in the same units than in mod_tides.F90 on hycom(cycles/day) */
  double freq[ ] =
  {
   0.893244048,0.929535695,0.997262079,1.002737898,1.895981946,1.932273593,2.000000000,2.005475795
  };


   //printf("Enter model dimensions nx and ny:\n"); 
   //scanf(" %d%d",&nx,&ny);
   //printf("dimensions nx=%d and ny=%d \n",nx,ny); 


   printf("Retrieve grid dim from regional.grid.b \n");
   flatlon=fopen("regional.grid.b","r");
   if (flatlon==NULL) {
      printf("fes2mod,fopen:Pb while opening regional.grid.b!!!\n");
      return -1;
   }
   if ( fgets(line, 100, flatlon) != NULL ) {
      printf(line);
      sscanf(line,"%d",&nx);
   } else {
      printf("fes2mod:Pb while parsing regional.grid.b for nx\n");
      return -1;
   }

   if ( fgets(line, 100, flatlon) != NULL ) {
      printf(line);
      sscanf(line,"%d",&ny);
   } else {
      printf("fes2mod:Pb while parsing regional.grid.b for ny\n");
      return -1;
   }
   fclose(flatlon);
   printf("nx = %d , ny = %d\n",nx, ny);


     
   float         lon4[nx*ny],lat4[nx*ny],dep4[nx*ny];
   double        mlon[nx][ny],mlat[nx][ny],depths[nx][ny];
   double        ampha[2*(nx+ny)][2*nconst],dep[nx*ny],lat[nx*ny],lon[nx*ny];
   // add some correction when nearest FES grid point is not wet
   double        dep_jp[2*(nx+ny)];
   int           Tinter[2*(nx+ny)];
   int           Tcorrect;
   int           max(int a,int b);
   int           min(int a,int b);
   double        w0;
   /// by Xiejp in 19/10/2010
   int           ti,offset, n2drec ;
   char* fespath;
   double huge;

   // One record length in .a - files - we need to seek to these positions
   n2drec=((nx*ny+4095)/4096)*4096;
   huge =1e30;
   Tcorrect=0;  
   for (i=0; i<2*(nx+ny); i++) Tinter[i]=0;

  /* Get hycom model grid info: depth  using regional files  */
   fdepths=fopen("regional.depth.a","rb");
   if ((fdepths !=NULL) ) {
      fseek(fdepths,0,SEEK_SET); // Direct access fortran Starts at 0 
      ti=fread(&dep4,sizeof(float),nx*ny,fdepths);
      if(ti !=nx*ny) {
        printf("fes2mod,fseek:Pb while accessing Depth file!!!\n");
        return -1;
      }
      fclose(fdepths);
   } else {
      printf("fes2mod,fopen:Pb while opening regional.grid.a!!!\n");
        return -1;
   }


  /* Get hycom model grid info: latlon  using regional files  */
   flatlon=fopen("regional.grid.a","rb");
   if ((flatlon !=NULL) ) {
      fseek(flatlon,0,SEEK_SET); // Direct access fortran Starts at 0 
      ti=fread(&lon4,sizeof(float),nx*ny,flatlon);
      if(ti !=nx*ny) {
        printf("fes2mod,fseek:Pb while accessing grid file!!!\n");
        return -1;
      }
      fseek(flatlon,n2drec*4,SEEK_SET); // Direct access fortran Starts n2drec*4(bytes)
      ti=fread(&lat4,sizeof(float),nx*ny,flatlon);
      if(ti !=nx*ny) {
        printf("fes2mod,fseek:Pb while accessing grid file!!!\n");
        return -1;
      }
      fclose(flatlon);
   } else {
      printf("fes2mod,fopen:Pb while opening regional.grid.a!!!\n");
        return -1;
   }
   printf("fes2mod: Done reading regional.grid/depth.a\n");
  k=0;
  for (j=0; j<ny;j++) {
     for (i=0; i<nx;i++) {
          //printf("swap before: %f \n",dep4[k]);
          swapByte(&dep4[k],sizeof(float));
          swapByte(&lat4[k],sizeof(float));
          swapByte(&lon4[k],sizeof(float));
          //printf("swap after: %f \n",dep4[k]);
          if (dep4[k] < 0.5 * huge ) {
             depths[i][j]=(double)dep4[k];
          } else {
             depths[i][j]=(double) 0.;
          }
          mlat[i][j]  =(double)lat4[k];
          mlon[i][j]  =(double)lon4[k];
          k=k+1; 
        }
   }

  /* Keep as an option */
  if (0==1) {
     /* Get hycom model grid info: depth+latlon!                      */
     sprintf(fdepth_name,"depths%dx%d.uf",nx,ny);
     fdepths=fopen(fdepth_name,"rb");
     flatlon=fopen("newpos.uf","rb");
     if ((fdepths !=NULL) & (flatlon !=NULL) )
     {
            fseek(fdepths,0,SEEK_SET);
            //ti=ftell(fdepths);
            //printf("size file in byte: %d \n",ti);
            fseek(fdepths,4,SEEK_SET); // Skip fortran 4-byte header 
            ti=fread(&dep,sizeof(double),nx*ny,fdepths);
            if(ti !=nx*ny) 
            {
              printf("fes2mod,fseek:Pb while accessing Depth file!!!");
              return -1;
            }
            fclose(fdepths);

      //************************************************************************//      
            fseek(flatlon,0,SEEK_SET);
            fseek(flatlon,4,SEEK_SET); // Skip fortran 4-byte header 
            ti=fread(&lat,sizeof(double),nx*ny,flatlon);
            if(ti !=nx*ny) 
            {
              printf("fes2mod,fseek:Pb while accessing lat data in newpos.uf file!!!");
              return -1;
            }
            
            //ti=ftell(flatlon);
            //printf("pos in file: %d \n",ti);
            ti=fread(&lon,sizeof(double),nx*ny,flatlon);
            if(ti !=nx*ny) 
            {
              printf("fes2mod,fseek:Pb while accessing lon data in newpos.uf file!!!");
              return -1;
            }
            //ti=ftell(flatlon);
            //printf("pos in file: %d \n",ti);
            fclose(flatlon);
     }
     else
     {
        printf("Enter the right latlon and depth model grid file names!!!!\n");
        return -1;
     }
     k=0;
     for (j=0; j<ny;j++)
     {
        for (i=0; i<nx;i++)
           {
             swapByte(&dep[k],sizeof(double));
             swapByte(&lat[k],sizeof(double));
             swapByte(&lon[k],sizeof(double));
             depths[i][j]=dep[k];
             mlat[i][j]=lat[k];
             mlon[i][j]=lon[k];
             k=k+1; 
           }
      }
      //printf("TOTOTO depth   %f %g %g\n",depths[249][223],mlat[98][224],mlon[180][117]);
   } /* if false */


 /*--------- Initialize memory for FES algorithms ------------------*/
 /* KAL - added option to specify FES_PATH */
  if ( getenv("FES_PATH") == NULL ) {
     printf("\nError: Set enironment variable FES_PATH to the location of FES data \n");
     return -1;
  } else {
     fespath=getenv("FES_PATH") ;
  }
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

 for (j=0;j<2*(nx+ny);j++) {
   for (i=0;i<2*nconst;i++) ampha[j][i]=-999.9;
 }

 for (j=0;j<nx;j++)
 {
    m=j;                                //South  0 <m<nx
    k=1;
    dlon=mlon[j][k];
    dlat=mlat[j][k];
    dep_jp[m]=depths[j][k];
    if (depths[j][k]>0)
    {
       rc=fesCore8Const(shortTide,dlat,dlon,nx,ny,nconst,m,&ampha);
       if ( rc != FES_SUCCESS ) 
       {
          Tinter[m]=1;
          Tcorrect=Tcorrect+1;
          printf("%s %f %f %f\n","     Error at South: ",dlon,dlat,depths[j][k]);
         // goto onError;  
       }
    }
    else
       for (i=0;i<2*nconst;i++) ampha[m][i]=-999.9;

    
    m=2*nx+ny-j-1;                   //North  jj+kk+1<m<2*jj+kk
    k=ny-2;
    dlon=mlon[j][k];
    dlat=mlat[j][k];
    dep_jp[m]=depths[j][k];
    if (depths[j][k]>0)
    {
       rc=fesCore8Const(shortTide,dlat,dlon,nx,ny,nconst,m,&ampha);
       if ( rc != FES_SUCCESS )
       {
          Tinter[m]=1;
          Tcorrect=Tcorrect+1;
          printf("%s %f %f %f\n","    Error at North: ",dlon,dlat,depths[j][k]);
         // goto onError;  
       }
    }
    else
       for (i=0;i<2*nconst;i++) ampha[m][i]=-999.9;
    //printf("2222ampha: %f %f %f %f\n",ampha[m][4],depths[j][k],dlat,dlon);
 }  

 for (k=0;k<ny;k++)
 {
    m=nx+k;                        //East     jj+1 <m<jj+kk
    j=nx-2;
    dlon=mlon[j][k];
    dlat=mlat[j][k];
    dep_jp[m]=depths[j][k];
    if (depths[j][k]>0)
    {
       rc=fesCore8Const(shortTide,dlat,dlon,nx,ny,nconst,m,&ampha);
       if ( rc != FES_SUCCESS )
       {
          Tinter[m]=1;
          Tcorrect=Tcorrect+1;
          printf("%s %f %f %f\n","     Error at East: ",dlon,dlat,depths[j][k]);
         // goto onError;  
       }
    }
    else
       for (i=0;i<2*nconst;i++) ampha[m][i]=-999.9;
   
    //printf("3333ampha: %f %f %f %f\n",ampha[m][4],depths[j][k],dlat,dlon);

    m=2*(nx+ny)-1-k;                   //West    2*jj+kk+1<m<2*(jj+kk)
    j=1;
    dlon=mlon[j][k];
    dlat=mlat[j][k];
    dep_jp[m]=depths[j][k];
    if (depths[j][k]>0)
    {
       rc=fesCore8Const(shortTide,dlat,dlon,nx,ny,nconst,m,&ampha);
       if ( rc != FES_SUCCESS )
       {
          Tinter[m]=1;
          Tcorrect=Tcorrect+1;
          printf("%s %f %f %f\n","     Error at West: ",dlon,dlat,depths[j][k]);
         // goto onError;  
       }
    }
    else
       for (i=0;i<2*nconst;i++) ampha[m][i]=-999.9;

    //printf("4444ampha: %f %f %f %f\n",ampha[m][4],depths[j][k],dlat,dlon);
    
 }

  // correct the elevation when model point is dry in FES
  if (Tcorrect>0) 
  {
     for (i=0;i<2*(nx+ny);i++) 
     {
        if (Tinter[i]==1)
        {
           printf("%s %d %s\n","     correct the grid ",i," !");
           j=max(i-1,0); 
           k=min(i+1,2*(nx+ny)); 
           //printf("%d %d %d\n",j,i,k);
           if (dep_jp[j]>0.&Tinter[j]==0)
           {
             if (dep_jp[k]>0.&Tinter[k]==0) 
               w0=0.5;
             else 
               w0=1.0;
           }
           else
           {
             if (dep_jp[k]>0.&Tinter[k]==0) 
               w0=0.;
             else 
               w0=1.0;
           }
           for (m=0;m<2*nconst;m++) 
              ampha[i][m]=w0*ampha[j][m]+(1.-w0)*ampha[k][m];
        }
     }
  }
  printf("%s\n",fespath);


/**************************   DUMP DATA IN FESobc_elev.dat  *********************/

 fout=fopen("FESobc_elev.dat","w+");

 fprintf(fout,"%s",var);
 for (i=0;i<nconst;i++)
 {
    fprintf(fout,",\"a%s\"",namewave[i]);
    fprintf(fout,",\"p%s\"",namewave[i]);
 }
 fprintf(fout,"\n");
 itec=96;
 fprintf(fout,"TEXT X=80, Y=%d, H=2,T=\"Model Dim.   %d %d\"\n",itec,nx,ny);
 itec=itec-3;
 fprintf(fout,"TEXT X=80, Y=%d, H=2,T=\"Constituents %d\"\n",itec,nconst);
 itec=itec-3;
 fprintf(fout,"TEXT X=80, Y=%d, H=2,T=\"Name Freq (cyc/d)\"\n",itec);
 for (i=0;i<nconst;i++)
 {
    k=iwave[i];
    itec=itec-3;
    fprintf(fout,"TEXT X=80, Y=%d, H=2,T='%s %10.9f'\n",itec,namewave[i],freq[i]);
 }
 fprintf(fout,"ZONE I= %d, F=POINT\n",2*(nx+ny));
 for (m=0;m<2*(nx+ny);m++)
 {
    fprintf(fout,"%4d",m+1);
    for (i=0;i<2*nconst;i++) fprintf(fout,"%8.2f",ampha[m][i]);
    fprintf(fout,"\n");
 }

 printf("\n");
 printf(" SW m=  1 \n"); 
 printf(" SE m=  %d\n",nx); 
 printf(" NE m=  %d\n",nx+ny); 
 printf(" NW m=  %d\n",2*(nx+ny)); 
 
 fprintf(fout,"GEOMETRY X=0, Y=0, T=LINE C=BLACK, CS=GRID\n");
 m=3;
 fprintf(fout,"%d\n",m);                // Number of lines
 m=2;                             
 i1=0;                             // On axis value
 i2=1000;                          // A max amp/phase value
 fprintf(fout,"%d\n",m);                // Points in line
 fprintf(fout,"%5d %5d\n",nx,i1);
 fprintf(fout,"%5d %5d\n",nx,i2);
 fprintf(fout,"%d\n",m);               // Points in line
 fprintf(fout,"%5d %5d\n",nx+ny,i1);
 fprintf(fout,"%5d %5d\n",nx+ny,i2);
 fprintf(fout,"%d\n",m);               // Points in line
 fprintf(fout,"%5d %5d\n",2*nx+ny,i1);
 fprintf(fout,"%5d %5d\n",2*nx+ny,i2);

 printf("\n");
 printf(" OBS data are dumped in FESobc_elev.dat\n");
 printf("\n");
 printf(" Note that FESobc_elev.dat is dumped in Tecplot format\n");


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
