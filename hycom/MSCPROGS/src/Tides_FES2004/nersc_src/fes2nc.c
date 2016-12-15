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
#include <netcdf.h>
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

void handle_error(int status) {
if (status != NC_NOERR) {
  fprintf(stderr, "%s\n", nc_strerror(status));
  exit(-1);
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
  int           nx,ny;
  int           i,j,k,m,itec,i1,i2; 
  double        dlon,dlat;
  FILE*         flatlon=NULL;
  FILE*         fdepths=NULL;
  char          fdepth_name[16];
  char          var[]="VARIABLES='m'"; 
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
   int           ti,offset, n2drec ;
   char* fespath;
   double huge;

   // One record length in .a - files - we need to seek to these positions
   n2drec=((nx*ny+4095)/4096)*4096;
   huge =1e30;

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
      printf("fes2mod,fopen:Pb while opening regional.depth.a!!!\n");
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
          if (dep4[k] < 0.5 * huge & dep4[k] > 0.) {
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
   j=520;  i=344;
   printf("TOTOTO depth   %f %g %g\n",depths[i][j],mlat[i][j],mlon[i][j]);



 /*--------- Initialize memory for FES algorithms ------------------*/
 /* KAL - added option to specify FES2004_PATH */
  if ( getenv("FES2004_PATH") == NULL ) {
     printf("\nError: Set enironment variable FES2004_PATH to the location of FES data \n");
     return -1;
  } else {
     fespath=getenv("FES2004_PATH") ;
  }
  //rc = fesNew(&shortTide, FES_TIDE, FES_IO, "/home/fimm/nersc/intissar/Progs/Tide_FES/data");
  rc = fesNew(&shortTide, FES_TIDE, FES_IO, fespath);
  if ( rc != FES_SUCCESS )
      goto onError;


/**************************   DUMP DATA IN netcdf file  *********************/
  int  ncid,dimids[2],varid, xdimid, ydimid, iconst,varamp[nconst] , varpha[nconst];
  static const char* netcdffile="fes2004.nc";
  char wname[3];
  char vname[100];
  size_t start[2] = {0, 0};
  size_t count[2] = {nx, ny};
  double    amphb[2*nconst];
  double amplitude[nx][ny][nconst]; 
  double phase[nx][ny][nconst]; 
  double tmpamp[nx][ny]; 
  double tmppha[nx][ny]; 
  static double fillval[]= {-999.9};
  // add some correction for bathymetry mismatch between and FES 
  int    ip0,ip1,jp0,jp1;
  int    Tinter[nx][ny];
  int    Ncorrect,Nout;
  int    max(int a, int b);
  int    min(int a, int b);
  double    tmp_a0[2*nconst];
  double    tmp_b0[2*nconst];
  int       tmp_n;
  Ncorrect=0;
  Nout=0;
  for (i=0;i<nx;i++) {
     for (j=0;j<ny;j++) {
        Tinter[i][j]=0;
     }
  }

  // by Xiejp in 27/10/2010
  

  handle_error(nc_create(netcdffile, NC_CLOBBER, &ncid));
  handle_error(nc_def_dim(ncid, "nx", (long) nx ,&xdimid ));
  handle_error(nc_def_dim(ncid, "ny", (long) ny ,&ydimid ));
  //dimids[1]=xdimid; dimids[0]=ydimid; 
  dimids[0]=xdimid; dimids[1]=ydimid; 

  handle_error(nc_def_var(ncid, "longitude" , NC_DOUBLE, 2, dimids , &varid) );
  handle_error(nc_enddef(ncid));
  handle_error(nc_put_vara_double(ncid, varid, start, count,  &mlon[0][0]  ));

  handle_error(nc_redef(ncid));
  handle_error(nc_def_var(ncid, "latitude" , NC_DOUBLE, 2, dimids , &varid) );
  handle_error(nc_enddef(ncid));
  handle_error(nc_put_vara_double(ncid, varid, start, count,  &mlat[0][0]  ));

  handle_error(nc_redef(ncid));
  handle_error(nc_def_var(ncid, "depth" , NC_DOUBLE, 2, dimids , &varid) );
  handle_error(nc_enddef(ncid));
  handle_error(nc_put_vara_double(ncid, varid, start, count,  &depths[0][0]  ));

  printf("Putting constituents in %s\n",netcdffile);
  printf("NB: somme errors usually occur ");
  for (iconst=0;iconst<nconst;iconst++) {
     handle_error(nc_redef(ncid));
     strcpy(wname,namewave[iconst]);
     wname[2]='\0';
     strcpy(vname,wname); strcat(vname,"_amplitude");
     handle_error(nc_def_var(ncid, vname, NC_DOUBLE, 2, dimids , &varamp[iconst]) );
     strcpy(vname,wname); strcat(vname,"_phase");
     handle_error(nc_def_var(ncid, vname, NC_DOUBLE, 2, dimids , &varpha[iconst]) );
     handle_error(nc_put_att_double (ncid, varamp[iconst], "_FillValue", NC_DOUBLE, 1, fillval));
     handle_error(nc_put_att_double (ncid, varpha[iconst], "_FillValue", NC_DOUBLE, 1, fillval));
     handle_error(nc_enddef(ncid));
  }

  /* Define constituent variables */
  for (i=0;i<nx;i++) {
     for (j=0;j<ny;j++) {
        //printf(" %d %d \n",i,j);
        if (depths[i][j]>0.1) {
           dlon=mlon[i][j];
           dlat=mlat[i][j];
           rc = diag8Const(shortTide,dlat,dlon,nconst,&amphb);
           if ( rc != FES_SUCCESS ) {
              //goto onError;
              for (iconst=0;iconst<2*nconst;iconst++) amphb[iconst]=-999.9;
              Tinter[i][j]=1;
              Ncorrect++;
              printf("#ERROR :%d  %d  %s\n", i,j, fesErr2String(rc));
           } 
           else {
              Tinter[i][j]=2;
           }
        } else {
           for (iconst=0;iconst<2*nconst;iconst++) amphb[iconst]=-999.9;
        }
        for (iconst=0;iconst<nconst;iconst++){
           amplitude[i][j][iconst]=amphb[2*iconst  ];
           phase    [i][j][iconst]=amphb[2*iconst+1];
        }
     }
  }

  // considering the correction by nearest value
  while (Ncorrect>0) {
    Nout++;
    for (i=0; i<nx; i++) {
      for (j=0; j<ny; j++) {
        if (Tinter[i][j]==1) {
          ip0=max(i-1,0);
          ip1=min(i+1,nx-1);
          jp0=max(j-1,0);
          jp1=min(j+1,ny-1);
          tmp_n=0;  
          for (iconst=0; iconst<2*nconst; iconst++) tmp_b0[iconst]=0;
          if (Tinter[ip0][jp0]==2) {
            tmp_n=tmp_n+1;
            for (iconst=0;iconst<nconst;iconst++){
              tmp_a0[2*iconst]=amplitude[ip0][jp0][iconst];
              tmp_a0[2*iconst+1]=phase[ip0][jp0][iconst];
            }
            for (iconst=0; iconst<2*nconst; iconst++) {
              tmp_b0[iconst]=(tmp_b0[iconst]*(tmp_n-1)+tmp_a0[iconst])/tmp_n;
            }
          //  printf(" %d %d %d %f %f\n",i,j,tmp_n,tmp_b0[1],tmp_a0[1]);
          }
          if (Tinter[ip1][jp1]==2) {
            tmp_n=tmp_n+1;
            for (iconst=0;iconst<nconst;iconst++){
              tmp_a0[2*iconst]=amplitude[ip1][jp1][iconst];
              tmp_a0[2*iconst+1]=phase[ip1][jp1][iconst];
            }
            for (iconst=0; iconst<2*nconst; iconst++) {
              tmp_b0[iconst]=(tmp_b0[iconst]*(tmp_n-1)+tmp_a0[iconst])/tmp_n;
            }
          //  printf(" %d %d %d %f %f\n",i,j,tmp_n,tmp_b0[1],tmp_a0[1]);
          }
          if (Tinter[ip0][jp1]==2) {
            tmp_n=tmp_n+1;
            for (iconst=0;iconst<nconst;iconst++){
              tmp_a0[2*iconst]=amplitude[ip0][jp1][iconst];
              tmp_a0[2*iconst+1]=phase[ip0][jp1][iconst];
            }
            for (iconst=0; iconst<2*nconst; iconst++) {
              tmp_b0[iconst]=(tmp_b0[iconst]*(tmp_n-1)+tmp_a0[iconst])/tmp_n;
            }
           // printf(" %d %d %d %f %f\n",i,j,tmp_n,tmp_b0[1],tmp_a0[1]);
          }
          if (Tinter[ip1][jp0]==2) {
            tmp_n=tmp_n+1;
            for (iconst=0;iconst<nconst;iconst++){
              tmp_a0[2*iconst]=amplitude[ip1][jp0][iconst];
              tmp_a0[2*iconst+1]=phase[ip1][jp0][iconst];
            }
            for (iconst=0; iconst<2*nconst; iconst++) {
              tmp_b0[iconst]=(tmp_b0[iconst]*(tmp_n-1)+tmp_a0[iconst])/tmp_n;
            }
          //  printf(" %d %d %d %f %f\n",i,j,tmp_n,tmp_b0[1],tmp_a0[1]);
          }
          // return the search value into the grid (ip,jp)
          if (tmp_n) {
            for (iconst=0;iconst<nconst;iconst++){
              amplitude[i][j][iconst]=tmp_b0[2*iconst  ];
              phase    [i][j][iconst]=tmp_b0[2*iconst+1];
            }
            Tinter[i][j]=2; 
          }

        }
      }
    }

    Ncorrect=0;
    for (i=0; i<nx; i++) {
      for (j=0; j<ny; j++) {
        if (Tinter[i][j]==1) Ncorrect++;
      }
    }
    printf("%s %d %s\n","Ncorrect=",Ncorrect," with iteration of ",Nout);
    // Breakout due to cannot correct it full
    if (Nout>20) { 
      Ncorrect=0;
      for (i=0; i<nx; i++) {
        for (j=0; j<ny; j++) {
          if (Tinter[i][j]==1) {
            printf("%s %d %d\n","Don't know: ",i,j);
          }
        }
      }
    }
  }


  for (iconst=0;iconst<nconst;iconst++) {
     printf("%s\n",namewave[iconst]);
     for (i=0;i<nx;i++) {
        for (j=0;j<ny;j++) {
           //printf("%d %d\n",i,j);
           tmpamp[i][j]=amplitude[i][j][iconst];
           tmppha[i][j]=phase    [i][j][iconst];
        }
     }
     handle_error(nc_put_vara_double(ncid, varamp[iconst], start, count,  &tmpamp[0][0]));
     handle_error(nc_put_vara_double(ncid, varpha[iconst], start, count,  &tmppha[0][0]));
  }
  handle_error(nc_close(ncid));
  printf("Constituents in %s\n",netcdffile);



  goto onTerminate;

onError:
  printf("#ERROR : %s\n", fesErr2String(rc)),
  rc = 1;

onTerminate:
  /* Free memory for FES algorithms */
  fesDelete(shortTide);
  return rc;
}
