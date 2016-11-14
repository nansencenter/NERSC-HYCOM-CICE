      SUBROUTINE lherov(uw,vw,mxlat,mxlon,mlat,mlon,
     +           iflg,nx,ny,mean)
c ---------------------------------------------------------
c --- Author: Knud Simonsen, NERSC,6/12 1993.
c ----------------------------------------------------------
c --- Reads the surface winds from Hellerman and Rosenstein
c --- and interpolate into scalar points.
c ----------------------------------------------------------
c --- External routines:
c --- rotate, fillup, intpol, rHeRo 
c ----------------------------------------------------------
c --- Input: mxlat,mxlon: The unit vector (0,1) in grid coordinates
c ---            in geogr. coordinates.
c ---        mlat,mlon: Geogr. coord. of the model grid points.
c ---            These input values are not used. New values 
c ---            are calc. by routine newold  .
c ---        iflg:  Land/sea mask
c ---        nx,ny: Size of the model domain.
c ---        mean:  mean or dummy   
c ----------------------------------------------------------------
      PARAMETER(nxd=180,nyd=90)
      INTEGER nx,ny,iflg(nx,ny)
      REAL uw(nx,ny,12),vw(nx,ny,12),mxlat(nx,ny),
     +     mxlon(nx,ny),mlat(nx,ny),mlon(nx,ny),a,b,
     +     du(nxd,nyd,12),dv(nxd,nyd,12),dlat(nyd),dlon(nxd)
      CHARACTER*15 filen
                                      !Read the wind field from
                                      !Hellerman and Rosenstein.
      CALL rHeRo(du,dv,dlat,dlon,nxd,nyd)
                                   !Interpolate U-comp
     
      CALL intpol(du,dlat,dlon,nxd,nyd,
     +        uw,mlat,mlon,nx,ny,12)

                                      !Interpolate V-comp
      CALL intpol(dv,dlat,dlon,nxd,nyd,
     +        vw,mlat,mlon,nx,ny,12)
                                      !Fillup missing values.
      CALL fillup(uw,nx,ny,iflg,12,mean)
      CALL fillup(vw,nx,ny,iflg,12,mean)
                                      !Rotate U-comp. V-comp
      CALL rotate(uw,vw,mxlat,mxlon,ny,nx)
c       write(*,*)'Not rotated'

!      write(*,*)'Perform additional interpolation at the Pole'
!      write(*,*)'around point 106,48'
!      DO i=103,109
!       DO j=45,51
!        DO m=1,12
!         uw(i,j,m)=-999.0
!         vw(i,j,m)=-999.0
!        ENDDO
!       ENDDO
!      ENDDO
!      CALL fillup(uw,nx,ny,iflg,12,mean)
!      CALL fillup(vw,nx,ny,iflg,12,mean)

                                      
      RETURN
      END

