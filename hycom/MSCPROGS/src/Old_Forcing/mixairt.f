      SUBROUTINE mixairt(ta,th,nx,ny,iflg,mlon,mlat)
c -------------------------------------------------------------
c --- Loads air temperature data from COADS. Origionally the
c --- COADS data does not contains data in ice  covered areas.
c --- In the present vesion data from Shea (1986) is used at
c --- high lat. (oberhuber, 1992). Here, the cloudcover data
c --- is used to blank out the areas origionally with no data
c --- in the COADS, and data from   ECMWF is included in
c --- central Arctic and coastal weather station data made
c --- available by the CDIAC are used to fill the gab between.
c --- This procedure was used by Simonsen and Haugan (JGR), 1996.
c -------------------------------------------------------------

      REAL ta(nx,ny,12),th(nx,ny,12),mlon(nx,ny),mlat(nx,ny)
     &     ,a,b,vmin,vmax

      INTEGER nx,ny, iflg(nx,ny)

      CHARACTER filen*15 

                                 !Load COADS air temp data and
                                 !interpolate them onto model grid
      CALL lcoads3(ta,th,nx,ny,mlon,mlat)       

c.ks      filen='airta.dat'                    !Dump to tec.
c.ks      CALL tec1(nx,ny,12,filen,mlon,mlat,ta,273.15)
      

      filen= 'Data/TSC'            !Read ECMWF temp
      a=10000./120.
      b= 273.16-75.
      vmin = 221.
      vmax = 321.
      CALL lecmw3(th,nx,ny,mlon,mlat,filen,a,b,vmin,vmax)

c.ks      filen='airte.dat'                    !Dump to tec.
c.ks      CALL tec1(nx,ny,12,filen,mlon,mlat,th,273.15)


      
      DO j=1,ny                    !Combine COADS and ECMWF data
       DO i=1,nx
        DO k=1,12
         fl = .5 + SIGN(.5,(ta(i,j,k)+900.))      !=1 if ta exist
c.diag         if(fl.LT..5)  write(*,*)ta(i,j,k),th(i,j,k),fl
         ta(i,j,k)= fl*ta(i,j,k)
     +               +(1.-fl)*th(i,j,k)
        ENDDO
       ENDDO
      ENDDO

c.ks      filen='airte.dat'                    !Dump to tec.
c.ks      CALL tec1(nx,ny,12,filen,mlon,mlat,ta,273.15)


                                   !Read data from obs. stations.
      CALL rstat(nx,ny,ta,mlon,mlat) 

c.ks      filen='airta1.dat'                    !Dump to tec.
c.ks      CALL tec1(nx,ny,12,filen,mlon,mlat,ta,273.15)
c.ks       STOP'test'

                                   !Fillup missing data
      CALL fillall(ta,nx,ny,12,271.0)
      RETURN
      END


 



