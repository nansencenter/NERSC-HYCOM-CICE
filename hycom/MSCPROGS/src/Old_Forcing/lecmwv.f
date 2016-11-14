      SUBROUTINE lecmwv(mu,mv,iflg,nxm,nym,mlon,mlat,
     +     mxlon,mxlat,filen,
     +     a,b,vmin,vmax,mean)
c -------------------------------------------------------------
c --- This routine reads data from from file in
c --- readcoads, interpolates them on the model grid
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
      PARAMETER(nxd=144,nyd=73)
      REAL mu(nxm,nym,12),mv(nxm,nym,12),
     +          mlat(nxm,nym),mlon(nxm,nym),
     +          du(nxd,nyd,12),dv(nxd,nyd,12),
     +          dlat(nyd),dlon(nxd),
     +          mxlon(nxm,nym),mxlat(nxm,nym)
      CHARACTER*15 filen
      REAL a,b,vmin,vmax,mean
      INTEGER nxm,nym,iflg(nxm,nym) 
c      write(*,*)'Ckeck lecmwv'
c      GOTO 100

      CALL recmwv(du,dv,dlat,dlon,nxd,nyd,
     +               filen,a,b,vmin,vmax)
                                           !Interpolate
      CALL intpol(du,dlat,dlon,nxd,nyd,
     +        mu,mlat,mlon,nxm,nym,12)
                                           !Fillup missing values
      CALL fillup(mu,nxm,nym,iflg,12,mean) 
                                           !Interpolate v
      CALL intpol(dv,dlat,dlon,nxd,nyd,
     +        mv,mlat,mlon,nxm,nym,12)
                                           !Fillup missing values
      CALL fillup(mv,nxm,nym,iflg,12,mean) 
                                           !Rotate.       
 100  CALL rotate(mu,mv,mlat,mlon,nxm,nym)
c      write(*,*)'Not rotated'

!      write(*,*)'Perform additional interpolation at the Pole'
!      write(*,*)'around point 106,48'
!      DO i=103,109
!       DO j=45,51
!        DO m=1,12
!         mu(i,j,m)=-999.0
!         mv(i,j,m)=-999.0
!        ENDDO
!       ENDDO
!      ENDDO
!      CALL fillup(mu,nxm,nym,iflg,12,mean)
!      CALL fillup(mv,nxm,nym,iflg,12,mean)

      RETURN
      END
