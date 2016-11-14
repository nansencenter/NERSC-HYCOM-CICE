      SUBROUTINE lecmw3(mdata,nxm,nym,mlon,mlat,filen,
     +     a,b,vmin,vmax)
c -------------------------------------------------------------
c --- This routine reads data from file 'filen'
c --- and interpolates on the model grid
c --- The format of the file is the T106 - grid, see OPYC-manual
c --- Author: Knud Simonsen, NERSC, 2/12 1993
c --- 
c --- Input: filen: Filename
c ---        a,b:   Constants used to convert to physival units.
c ---        vmin,vmax: Range off accepted values 
c ---        nxm,nym: Size of the model grid.
c ---        mlat:  Contains lattitude of the model gridcells
c ---        mlon:  Contains longitude of the model gridcells. 
c ---        nyd:   Number of lattitudes. 
c ---        mean: A dummy  which is set forthe missing value
c ---              before the interpolation start
c --- Output: mdata: Datafile 
c -------------------------------------------------------------
      PARAMETER(nxd=320,nyd=160)
      REAL mdata(nxm,nym,12),mlat(nxm,nym),mlon(nxm,nym),
     +          dd(nxd,nyd,12),dlat(nyd),dlon(nxd)
      CHARACTER*15 filen
      REAL a,b,vmin,vmax,fl
      INTEGER nxm,nym,ok

      CALL recmw2(dd,dlat,dlon,nxd,nyd,
     +               filen,a,b,vmin,vmax)
      DO j=1,nyd                             !Pick Out data to be 
       DO i=1,nxd                            !applied in the Arctic
         ok = 0
         IF((j.GE.nyd-3).AND.(j.LE.nyd-2))ok=1 

         IF((j.GE.nyd-6.AND.j.LE.nyd-4).
     +      AND.(i.LE.27.OR.i.GE.257))ok=1


         IF((j.GE.nyd-8.AND.j.LE.nyd-7).
     +      AND.(i.LE.25.OR.i.GE.267))ok=1

         IF(j.GE.nyd-10.AND.j.LE.nyd-9.
     +      AND.(i.LE.24.OR.i.GE.277))ok=1

         IF(j.GE.nyd-12.AND.j.LE.nyd-11.
     +      AND.(i.LE.19.OR.i.GE.304))ok=1
         fl=FLOAT(ok)
         DO k=1,12
          dd(i,j,k)=fl*dd(i,j,k)- (1.-fl)*999.
         ENDDO 
       ENDDO
      ENDDO 
                                           !Interpolate
      CALL intpol(dd,dlat,dlon,nxd,nyd,
     +        mdata,mlat,mlon,nxm,nym,12,
     +        gridn,ypivn,xpivn,ypivo)

      RETURN
      END
