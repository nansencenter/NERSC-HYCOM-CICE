      SUBROUTINE lwalsh(mdata,iflg,nxm,nym,mlon,mlat,filen,
     +                  a,b,vmin,vmax)
c -------------------------------------------------------------
c --- This routine reads data from from file in
c --- readcoads, interpolates them on the model grid
c --- Author: Knud Simonsen, NERSC, 2/12 1993
c ---    modefied    24/2 1997
c --- 
c --- Input: filen: Filename
c ---        a,b:   Constants used to convert to physival units.
c ---        vmin,vmax: Range off accepted values 
c ---        nxm,nym: Size of the model grid.
c ---        mlat:  Contains lattitude of the model gridcells
c ---        mlon:  Contains longitude of the model gridcells. 
c ---        iflg:  Land/sea mask.
c ---        nyd:   Number of lattitudes. 
c ---        mean: A dummy  which is set forthe missing value
c ---              before the interpolation start
c --- Output: mdata: Datafile 
c -------------------------------------------------------------
      PARAMETER(nxd=80,nyd=58)
      REAL mdata(nxm,nym,12),mlat(nxm,nym),mlon(nxm,nym),
     +          dd(nxd,nyd,12),dlat(nxd,nyd),dlon(nxd,nyd)
      CHARACTER*20 filen
      REAL a,b,vmin,vmax
      INTEGER nxm,nym,iflg(nxm,nym) 

      write(*,*) filen
      CALL rwalsh(dd,dlat,dlon,nxd,nyd,
     +               filen,a,b,vmin,vmax)
      write(*,*) filen

                                           !Interpolate
      write(*,*) 'intpol2'
      CALL intpol2(dd,dlat,dlon,nxd,nyd,
     +        mdata,mlat,mlon,nxm,nym,12)
      write(*,*) 'intpol2'


      RETURN
      END
