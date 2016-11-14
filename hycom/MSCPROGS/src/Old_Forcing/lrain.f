      SUBROUTINE lrain(mdata,iflg,nxm,nym,mlon,mlat,filen,
     +     a,b,vmin,vmax)
c -------------------------------------------------------------
c --- This routine reads data from from file in
c --- readcoads, interpolates them on the model grid
c --- Author: Knud Simonsen, NERSC, 2/12 1993
c --- 
c --- Input: filen: Filename
c ---        a,b:   Constants used to convert to physival units.
c ---        vmin,vmax: Range off accepted values 
c ---        gridn:  Gridresolution
c ---        ypivo:  Model equator meridiane
c ---        xpivn:  i-index for the ME
c ---        ypivn:  j-index for true eq.
c ---        nxm,nym: Size of the model grid.
c ---        mlat:  Contains lattitude of the model gridcells
c ---        mlon:  Contains longitude of the model gridcells. 
c ---        iflg:  Land/sea mask.
c ---        nyd:   Number of lattitudes. 
c --- Output: mdata: Datafile 
c -------------------------------------------------------------
      PARAMETER(nxd=360,nyd=180)
      REAL mdata(nxm,nym,12),mlat(nxm,nym),mlon(nxm,nym),
     +          dd(nxd,nyd,12),dlat(nyd),dlon(nxd)
      CHARACTER*15 filen
      REAL a,b,vmin,vmax,gridn,ypivn,xpivn,ypivo
      INTEGER nxm,nym,iflg(nxm,nym) 

      CALL rrain(dd,dlat,dlon,nxd,nyd,
     +               filen,a,b,vmin,vmax)
c       WRITE(2,*)'TITLE=""'
c       WRITE(2,*)'VARIABLES="lon","LAT","P"'
c       WRITE(2,*)'ZONE I=360,J=180,K=12,F=BLOCK'
c       WRITE(2,'(10f10.3)')(((dlon(i),i=1,nxd),j=1,nyd),k=1,12)
c       WRITE(2,'(10f10.3)')(((dlat(j),i=1,nxd),j=1,nyd),k=1,12)
c       WRITE(2,'(10f10.3)')(((dd(i,j,k),i=1,nxd),j=1,nyd),k=1,12)
c       STOP
   
                                           !Interpolate
      CALL intpol(dd,dlat,dlon,nxd,nyd,
     +        mdata,mlat,mlon,nxm,nym,12)
      RETURN
      END
