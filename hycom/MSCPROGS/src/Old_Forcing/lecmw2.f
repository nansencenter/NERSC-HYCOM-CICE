      SUBROUTINE lecmw2(mdata,iflg,nxm,nym,mlon,mlat,filen,
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
c ---        iflg:  Land/sea mask.
c ---        nyd:   Number of lattitudes. 
c ---        mean: A dummy  which is set forthe missing value
c ---              before the interpolation start
c --- Output: mdata: Datafile 
c -------------------------------------------------------------
      PARAMETER(nxd=320,nyd=160)
      REAL mdata(nxm,nym,12),mlat(nxm,nym),mlon(nxm,nym),
     +          dd(nxd,nyd,12),dlat(nyd),dlon(nxd)
      CHARACTER*15 filen
      REAL a,b,vmin,vmax
      INTEGER nxm,nym,iflg(nxm,nym) 

      CALL recmw2(dd,dlat,dlon,nxd,nyd,
     +               filen,a,b,vmin,vmax)
                                           !Interpolate
      CALL intpol(dd,dlat,dlon,nxd,nyd,
     +        mdata,mlat,mlon,nxm,nym,12)

      RETURN
      END
