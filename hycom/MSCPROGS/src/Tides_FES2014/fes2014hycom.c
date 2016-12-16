/*
 *  C program for the FES prediction software. Adapted to read hycom grid, and provide
 *  tidal predictions for times given as input to this routine. times are relative
 *  to 1950-01-01 00:00:00 UTC
 *
 *  File      : fes2hycnc.c
 *  Date      : December 2016 
 *  Added by KAL
 *
 *  This prog reads info about the model grid lon-lat-depth-dimensions. Then calculates
 *  tidal time series for the model grid and writes the output to a netcdf file.
 *
 *  It reads the following files:
 *   regional.grid.a   --- contains model lon lat positions. (NB: big-endian)
 *   regional.depth.a  --- contains model bathymetry.        (NB: big-endian)
 *   
 * ----------------------------------------------------------------------
 *
 *  time  Julian day (days since 1950-01-01 00:00:00.000 UTC).
 *
 */       

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <string.h>

#include "fes.h"

/*
// ///////////////////////////////////////////////////////////////////////////
// function swaps a "size" byte variable, i.e. it converts their
// representation between little and big endian (in either direction).
//
// Parameters:
//   p:		pointer to the variable to swap
//   n:		size of the variable
// ///////////////////////////////////////////////////////////////////////////
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

// Writes error string from netcdf status and exits
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
   Main program
*/
int main(int argc, char **argv)
{ 
   FES shortTide;
  int		rc	= 0;
  int           nx,ny;
  int           i,j,k;
  double        dlon,dlat;
  FILE*         flatlon=NULL;
  FILE*         fdepths=NULL;
  char          line[100];

  int  ncid,dimids[2],xdimid, ydimid;
  static const char* netcdffileu="fes2014ts_u.nc";
  static const char* netcdffilev="fes2014ts_v.nc";
  static const char* netcdffileh="fes2014ts_h.nc";
  static const char* netcdffile;
  size_t start[2]  ;//= {0, 0};
  size_t count[2] ;// = {ny, nx};
  static double fillval[]= {-999.9};
  int    im1,jm1, ip1,jp1;
  int    max(int a, int b);
  int    min(int a, int b);
  int    Nmissing, it;
  int  dimids3[3],varid_lon, varid_lat, varid_depth, varid_elev, tdimid, varid_time, dimidst[1];
  static char timeunit[] = "seconds since 1950-01-01 00:00:00";
  size_t start3[3];
  size_t count3[3];
  double       h,hlp;
  double  jday ;
  double myfillv=fillval[0];
  
  int           ti,n2drec ;
  char* fespath;
  double huge;
  int i2,j2, itime;
  size_t countt[1], startt[1];

  int ntimes=argc-2;
  double        tidetimes[ntimes];
  int make_central_island=1;
  int stripe;
  static const char * inifnameh = "/ocean_tide.ini";
  static const char * inifnameu = "/eastward_velocity.ini";
  static const char * inifnamev = "/northward_velocity.ini";
  static const char * inifname;
  static const char * varnameu = "u";
  static const char * varnamev = "v";
  static const char * varnameh = "h";
  static const char * varname;
  static const char * varunitu = "m s**-1";
  static const char * varunitv = "m s**-1";
  static const char * varunith = "m";
  static const char * varunit;
  static const char * varsnu = "eastward_sea_water_velocity";  // Well, not exactly, but cant find a suitable CF name
  static const char * varsnv = "northward_sea_water_velocity"; // Well, not exactly, but cant find a suitable CF name
  static const char * varsnh = "sea_surface_height";           // Well, not exactly, but cant find a suitable CF name
  static const char * varsn;

  if (argc < 3) {
     printf("Incorrect arguments tide times provided\n");
     printf("Usage for tidal elevations: fes2hycom h time1 time2 time3 .... \n");
     printf("Usage for tidal current(u): fes2hycom u time1 time2 time3 .... \n");
     printf("Usage for tidal current(v): fes2hycom v time1 time2 time3 .... \n");
     exit(1);
  }

  printf("%d\n",strcmp(argv[1],"h"));
  if (! strcmp(argv[1],"h") )
  {
     printf("Computing tidal elevation(h)\n");
     netcdffile=netcdffileh;
     inifname=inifnameh;
     varname=varnameh;
     varunit=varunith;
     varsn=varsnh;
  }
  else if (! strcmp(argv[1],"u") )
  {
     printf("Computing tidal current eastwards(u)\n");
     netcdffile=netcdffileu;
     inifname=inifnameu;
     varname=varnameu;
     varunit=varunitu;
     varsn=varsnu;
  }
  else if (!strcmp(argv[1],"v") )
  {
     printf("Computing tidal current northwards(v)\n");
     netcdffile=netcdffilev;
     inifname=inifnamev;
     varname=varnamev;
     varunit=varunitv;
     varsn=varsnv;
  }
  else
  {
     printf("dont know what tidal component to calculate. \n");
     printf("Usage for tidal elevations: fes2hycom h time1 time2 time3 .... \n");
     printf("Usage for tidal current(u): fes2hycom u time1 time2 time3 .... \n");
     printf("Usage for tidal current(v): fes2hycom v time1 time2 time3 .... \n");
     exit(0);
  }



  


  //  Read command line args
  for (i = 2; i < argc; ++i)
  {
      tidetimes[i-1]=strtod(argv[i],NULL);

      // if first argument is negative, calculate tides for full domain
      if (i==1 && tidetimes[i-1] < 0.) 
      {
         make_central_island=0;
      }
      tidetimes[i-1]=fabs( tidetimes[i-1]);
      //printf("%f\n",tidetimes[i-1]);
  }


   // Get hycom grid dimensions from regional.grid.b
   printf("Retrieve grid dim from regional.grid.b \n");
   flatlon=fopen("regional.grid.b","r");
   if (flatlon==NULL) {
      printf("fes2hycnc,fopen:Pb while opening regional.grid.b!!!\n");
      return -1;
   }
   if ( fgets(line, 100, flatlon) != NULL ) {
      printf(line);
      sscanf(line,"%d",&nx);
   } else {
      printf("fes2hycnc:Pb while parsing regional.grid.b for nx\n");
      return -1;
   }

   if ( fgets(line, 100, flatlon) != NULL ) {
      printf(line);
      sscanf(line,"%d",&ny);
   } else {
      printf("fes2hycnc:Pb while parsing regional.grid.b for ny\n");
      return -1;
   }
   fclose(flatlon);
   printf("nx = %d , ny = %d\n",nx, ny);

   // KAL: Declare on stack. Problematic due to stack size issues.
   /*
   //float         lon4[nx*ny],lat4[nx*ny],dep4[nx*ny];
   //double        mlon[ny][nx],mlat[ny][nx],depths[ny][nx];
   //double        dep[nx*ny],lat[nx*ny],lon[nx*ny];
   //int           Tinter[ny][nx];
   //int           Tinteri[ny][nx];
   //int           Tinterj[ny][nx];
   */

   // malloc on heap. 
   float  *dep4     = malloc(nx*ny * sizeof(float));
   float  *lon4     = malloc(nx*ny * sizeof(float));
   float  *lat4     = malloc(nx*ny * sizeof(float));
   // This construct declares a contigous memory region for 2D arrays, but may
   // not be allowed on old compilers
   double (*mlon)   [nx] = malloc( sizeof(*mlon   )* ny);
   double (*mlat)   [nx] = malloc( sizeof(*mlat   )* ny);
   double (*depths) [nx] = malloc( sizeof(*depths )* ny);
   double (*tmpfld) [nx] = malloc( sizeof(*tmpfld )* ny);
   int    (*Tinter) [nx] = malloc( sizeof(*Tinter )* ny);
   int    (*Tinteri)[nx] = malloc( sizeof(*Tinteri)* ny);
   int    (*Tinterj)[nx] = malloc( sizeof(*Tinterj)* ny);


   // One record length in .a - files - we need to seek to these positions
   n2drec=((nx*ny+4095)/4096)*4096;
   huge =1e30;

   /* Get hycom model grid info: depth  using regional files  */
   fdepths=fopen("regional.depth.a","rb");
   if ((fdepths !=NULL) ) {
      fseek(fdepths,0,SEEK_SET); // Direct access fortran Starts at 0 
      ti=fread(dep4,sizeof(float),nx*ny,fdepths);
      if(ti !=nx*ny) {
        printf("fes2hycnc,fseek:Pb while accessing Depth file!!!\n");
        return -1;
      }
      fclose(fdepths);
   } else {
      printf("fes2hycnc,fopen:Pb while opening regional.depth.a!!!\n");
        return -1;
   }

   /* Get hycom model grid info: latlon  using regional files  */
   flatlon=fopen("regional.grid.a","rb");
   if ((flatlon !=NULL) ) {
      fseek(flatlon,0,SEEK_SET); // Direct access fortran Starts at 0 
      ti=fread(lon4,sizeof(float),nx*ny,flatlon);
      if(ti !=nx*ny) {
        printf("fes2hycnc,fseek:Pb while accessing grid file!!!\n");
        return -1;
      }
      fseek(flatlon,n2drec*sizeof(float),SEEK_SET); // Direct access fortran Starts n2drec*4(bytes)
      ti=fread(lat4,sizeof(float),nx*ny,flatlon);
      if(ti !=nx*ny) {
        printf("fes2hycnc,fseek:Pb while accessing grid file!!!\n");
        return -1;
      }
      fclose(flatlon);
   } else {
      printf("fes2hycnc,fopen:Pb while opening regional.grid.a!!!\n");
        return -1;
   }
   printf("fes2hycnc: Done reading regional.grid/depth.a\n");
  k=0;
  for (j=0; j<ny;j++) {
     for (i=0; i<nx;i++) {
          //printf("swap before: %f \n",dep4[k]);
          swapByte(&dep4[k],sizeof(float));
          swapByte(&lat4[k],sizeof(float));
          swapByte(&lon4[k],sizeof(float));
          //printf("swap after: %f \n",dep4[k]);
          if ((dep4[k] < 0.5 * huge) & (dep4[k] > 0.)) {
             depths[j][i]=(double)dep4[k];
          } else {
             depths[j][i]=(double) 0.;
          }
          mlat[j][i]  =(double)lat4[k];
          mlon[j][i]  =(double)lon4[k];
          //mlon[j][i]  =fmod(mlon[j][i]+360.,360.);
          //printf("%d %d %f\n",j,i,depths[j][i]);

          k=k+1; 
        }
   }
   j=10;  i=nx/2;
   printf("Test depth i=%d d=%d   %f %g %g\n",i,j,depths[j][i],mlat[j][i],mlon[j][i]);

  // Create an island inside domain. This way only a stripe along the boundary will contain tide data. 
  // //Remove if you need 
  // full output, but computations will take some time....
  stripe=50;
  if (make_central_island==1) 
  {
     printf("Warning: setting depths=0. in center of domain\n");
     for (j=0; j<ny;j++) 
     {
        for (i=0; i<nx;i++) 
        {
           if (i>stripe && i < nx-stripe && j>stripe && j < ny-stripe)
           {
              depths[j][i] = 0.;
           }
        }
     }
   }
      

  //--------- Initialize memory for FES algorithms ------------------
  if ( getenv("FES2014_PATH") == NULL ) {
     printf("Error: Set enironment variable FES2014_PATH to the location of FES data \n");
     return -1;
  } else {
     fespath=getenv("FES2014_PATH") ;
  }
  printf("FES2014_PATH set to %s\n",fespath);
  int newSize = strlen(fespath)  + strlen(inifname) + 1; 
  char * inifile = (char *)malloc(newSize);
  strcpy(inifile,fespath);
  strcat(inifile,inifname); // or strncat 
  printf("FES ini file is set to %s\n",inifile);

  // Open fes dataset
  //setenv("FES_DATA","/work/shared/nersc/msc/ModelInput/tides/FES2014/",1) ;
  //rc = fes_new(&shortTide, FES_TIDE, FES_IO,"/work/shared/nersc/msc/ModelInput/tides/FES2014/ocean_tide.ini");
  //size_t fes_access_mode=FES_IO;
  size_t fes_access_mode=FES_MEM;
  setenv("FES_DATA",fespath,1);
  if (fes_new(&shortTide, FES_TIDE, fes_access_mode,inifile)) 
  {
     printf("fes_new error : %s\n", fes_error(shortTide));
     fes_delete(shortTide);
     return 1;
  }
  else
  {
     printf("fes_new call returned success-ready to compute tides\n");
  }

  // Set buffer size to 1024 MB if FES_IO is chosen
  if (fes_access_mode == FES_IO) {
     if (fes_set_buffer_size(shortTide,1024))
     {
        printf("fes_set_buffer_size error : %s\n", fes_error(shortTide));
        fes_delete(shortTide);
        return 1;
     }
  }



  // *************************   DUMP DATA IN netcdf file  *********************


  // Create netcf file wfile with time series
  start[0]=0;
  start[1]=0;
  count[0] = ny;
  count[1] = nx;
  handle_error(nc_create(netcdffile, NC_CLOBBER, &ncid));
  handle_error(nc_def_dim(ncid, "nx", (long) nx ,&xdimid ));
  handle_error(nc_def_dim(ncid, "ny", (long) ny ,&ydimid ));
  handle_error(nc_def_dim(ncid, "time", NC_UNLIMITED ,&tdimid ));
  dimids[1]=xdimid; dimids[0]=ydimid; 
  dimids3[0]=tdimid; dimids3[1]=ydimid; dimids3[2]=xdimid; 
  dimidst[0]=tdimid;
  //
  handle_error(nc_def_var(ncid, "longitude" , NC_DOUBLE, 2, dimids  , &varid_lon) );
  handle_error(nc_def_var(ncid, "latitude"  , NC_DOUBLE, 2, dimids  , &varid_lat) );
  handle_error(nc_def_var(ncid, "depth"     , NC_DOUBLE, 2, dimids  , &varid_depth) );
  handle_error(nc_def_var(ncid, "time"      , NC_DOUBLE, 1, dimidst  , &varid_time) );
  handle_error(nc_def_var(ncid, varname      , NC_DOUBLE, 3, dimids3 , &varid_elev) );
  handle_error(nc_put_att_double (ncid, varid_elev, "_FillValue", NC_DOUBLE, 1, fillval));
  handle_error(nc_put_att_text   (ncid, varid_elev, "units", strlen(varunit),  varunit));
  handle_error(nc_put_att_text   (ncid, varid_elev, "standard_name", strlen(varsn),  varsn));
  handle_error(nc_put_att_text   (ncid, varid_time, "units", strlen(timeunit), timeunit));
  handle_error(nc_enddef(ncid));
  printf("Putting timeseries in %s\n",netcdffile);
  //
  handle_error(nc_put_vara_double(ncid, varid_lon, start, count,   &mlon[0][0]  ));
  handle_error(nc_put_vara_double(ncid, varid_lat, start, count,   &mlat[0][0]  ));
  handle_error(nc_put_vara_double(ncid, varid_depth, start, count, &depths[0][0]  ));




  // First pass. Get Time series and set Tinter flag
  printf("Getting ocean points where FES returns data\n");
  Nmissing=0;
  jday=tidetimes[0];
  for (i=0;i<nx;i++) {
     for (j=0;j<ny;j++) {
        Tinter[j][i]=0;
        if (depths[j][i]>0.1) {
           dlon=mlon[j][i];
           dlat=mlat[j][i];
           if (fes_core(shortTide,dlat,dlon,jday,&h,&hlp) ) 

           {
              printf("#FES ERROR  at i=%5d  j=%5d  lat=%10f lon=%10f %s\n", i,j,dlat,dlon, fes_error(shortTide)); 
              Tinter[j][i]=1;
           } 
           else
           {
              Tinter[j][i]=0;
              Tinteri[j][i]=i;
              Tinterj[j][i]=j;
              //printf(" %d %d %f\n",i,j,h);
           }
        }
        if (Tinter[j][i] == 1) 
        {
           Nmissing=Nmissing+1;
        }
     }
  }

  // Use nearest neighbour for missing points
  if (Nmissing > 0 )  
  {
     printf("Getting points where FES needs to be extrapolated to fit bathymetry\n");
  }
  else
  {
     printf("FES was defined for all points we wanted. No need for extraploation :-)\n");
  }
  it=1;
  while (Nmissing > 0 && it < 150) 
  {
     printf("Iteration %d, Nmissing=%d\n",it,Nmissing);
     Nmissing=0;
     for (i=0;i<nx;i++) {
        for (j=0;j<ny;j++) {
           if (depths[j][i]>0.1 && Tinter[j][i]==1) {
              ip1=min(nx-1,i+1);
              jp1=min(ny-1,j+1);
              im1=max(0,i-1);
              jm1=max(0,j-1);

              if (depths[j][ip1] > .1 && Tinter[j][ip1]==0) 
              {
                 Tinter[j][i]=0;
                 Tinteri[j][i]=Tinteri[j][ip1];
                 Tinterj[j][i]=Tinterj[j][ip1];
                 printf(" (%d,%d) -> (%d,%d) \n",ip1,j,Tinteri[j][ip1],Tinterj[j][ip1]);
              }
              else if (depths[j][im1] > .1 && Tinter[j][im1]==0) 
              {
                 Tinter[j][i]=0;
                 Tinteri[j][i]=Tinteri[j][im1];
                 Tinterj[j][i]=Tinterj[j][im1];
                 printf(" (%d,%d) -> (%d,%d) \n",im1,j,Tinteri[j][im1],Tinterj[j][im1]);
              }
              else if (depths[jm1][i] > .1 && Tinter[jm1][i]==0) 
              {
                 Tinter[j][i]=0;
                 Tinteri[j][i]=Tinteri[jm1][i];
                 Tinterj[j][i]=Tinterj[jm1][i];
                 printf(" (%d,%d) -> (%d,%d) \n",i,jm1,Tinteri[jm1][i],Tinterj[jm1][i]);
              }
              else if (depths[jp1][i] > .1 && Tinter[jp1][i]==0) 
              {
                 Tinter[j][i]=0;
                 Tinteri[j][i]=Tinteri[jp1][i];
                 Tinterj[j][i]=Tinterj[jp1][i];
                 printf(" (%d,%d) -> (%d,%d) \n",i,jp1,Tinteri[jp1][i],Tinterj[jp1][i]);
              }
           }
           if (Tinter[j][i] == 1) 
           {
              Nmissing+=1;
           }
        }
     }
     it+=1;
  }
  if (Nmissing > 0 )  
  {
     printf("Unable to extrapolate  FES  points\n");
     //fes_delete(shortTide);
     //return 1;
  }
  else
  {
     printf("Extrapolation successful\n");
  }

  // Calculate remaini pass. Get Time series and set Tinter flag
  for (itime=0;itime<ntimes;itime++) {
     jday=tidetimes[itime];
     printf("time = %.2f\n",jday);
     for (i=0;i<nx;i++) {
        //printf("i = %d\n",i);
        for (j=0;j<ny;j++) {
           tmpfld[j][i]=myfillv;
           if (depths[j][i]>0.1) {
              i2=Tinteri[j][i];
              j2=Tinterj[j][i];
              dlon=mlon[j2][i2];
              dlat=mlat[j2][i2];
              //printf("%d %d %f %f\n",i2,j2,dlon,dlat);
              if (fes_core(shortTide,dlat,dlon,jday,&h,&hlp) ) 
              {
                 printf("#ERROR :%d  %d  %s\n", i,j, fes_error(shortTide)); 
              } 
              else
              {
                 tmpfld[j][i]=h*.01; // Conversion cm ->m 
              }
           }
        }
     }
     start3[0]=itime; start3[1]=0; start3[2]=0; 
     count3[0]=1; count3[1]=ny; count3[2]=nx; 
     startt[0]=itime;
     countt[0]=1; 
     jday=jday*86400.;
     handle_error(nc_put_vara_double(ncid, varid_time, startt, countt,  &jday          ));
     handle_error(nc_put_vara_double(ncid, varid_elev, start3, count3,  &tmpfld[0][0]  ));
  }

  handle_error(nc_close(ncid));
  printf("Time series in %s\n",netcdffile);
  printf("Warning: Only a %d wide strip along the domain has tides. This is ok for barotropic nesting, but may not be what you need\n",stripe);
  printf("Warning: To actually calculate tides for the entire domain, make sure first time is negative.. NB: This can be slow \n");

  // Free memory for FES algorithms 
  fes_delete(shortTide);
  return rc;
}
