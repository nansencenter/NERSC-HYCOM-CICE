      SUBROUTINE lcoads3(mdata,th,nxm,nym,mlon,mlat)
c -------------------------------------------------------------
c --- This routine reads air temp data from COADS
c --- and use the CC field to blank the area in the
c --- airtemp-field, which oriogionally did not contained data.
c --- In addition not using COADS north of 63N except in the 
c --- Nordic Seas
c -------------------------------------------------------------
c --- Author: Knud Simonsen, NERSC, 2/12 1993
c --- 
c --- Input: filen: Filename
c ---        nxm,nym: Size of the model grid.
c ---        mlat:  Contains lattitude of the model gridcells
c ---        mlon:  Contains longitude of the model gridcells. 
c --- Output: mdata: Datafile 
c -------------------------------------------------------------
      PARAMETER(nxd=180,nyd=90)
      REAL mdata(nxm,nym,12),mlat(nxm,nym),mlon(nxm,nym),
     +          dd(nxd,nyd,12),th(nxd,nyd,12),
     +          dlat(nyd),dlon(nxd)
      CHARACTER*15 filen
      REAL a,b,vmin,vmax,gridn,ypivn,xpivn,ypivo,fl,
     +      fl1,fl2,fl3
      INTEGER nxm,nym,iw,ie,js,jn 
     
      filen= 'Data/AIRGLOB'    !Read COADS air temp. In the       
      a=100.                      !Arctic the data from Shea(1986)
      b = 273.15-60.              !and the Antarctic the data
      vmin = 221.
      vmax = 321.
      CALL rcoads(dd,dlat,dlon,nxd,nyd,filen,a,b,vmin,vmax)

      filen= 'Data/COVER'      !Read Cloud cover  
      a= 80.
      b= 0.
      vmin=0.
      vmax=.98
      CALL rcoads(th,dlat,dlon,nxd,nyd,filen,a,b,vmin,vmax)
c.diag      filen='airth.dat'                    !Dump to tec.
c.diag      CALL tecij(nxd,nyd,12,filen,th,0.)

      DO i=1,nxd                           !Perform blanking
       DO j=1,nyd
         fl1 = 1. 
         DO k=1,12
          fl = .5 + SIGN(.5,th(i,j,k))     !=1 if COADS exist
          fl1 = MIN(fl,fl1)
         ENDDO
c         th(i,j,1)=fl1                    !=1 if COADS data exsist
                                           !in all 12 months else = 0 
         DO k=1,12
           dd(i,j,k)= fl1*dd(i,j,k)- (1.-fl1)*999. 
         ENDDO
       ENDDO 
      ENDDO 
c.diag      write(*,'(20f6.1)')((dd(i,j,1),i=1,nxd),j=1,nyd)
c.diag      write(*,'(20f6.1)')((th(i,j,1),i=1,nxd),j=1,nyd)
        
                                        !In addition:
      DO j=73,nyd                       !Not using COADS north of 63
       DO i=1,nxd                       !except in the Nordic Seas
        DO k=1,12
         fl=dd(i,j,k)
         dd(i,j,k)=-999. 
         IF(j.LE.78.AND.i.GE.87.AND.i.LE.97)dd(i,j,k)=fl
         IF(j.EQ.79.AND.i.GE.87.AND.i.LE.97)dd(i,j,k)=fl
         IF(j.EQ.80.AND.i.GE.87.AND.i.LE.97)dd(i,j,k)=fl
         IF(j.EQ.81.AND.i.GE.88.AND.i.LE.97)dd(i,j,k)=fl
         IF(j.EQ.82.AND.i.GE.88.AND.i.LE.97)dd(i,j,k)=fl
c.ks         IF(j.EQ.83.AND.i.GE.88.AND.i.LE.96)dd(i,j,k)=fl
        ENDDO
       ENDDO
      ENDDO
                                        !Interpolate onto model grid
      CALL intpol(dd,dlat,dlon,nxd,nyd,
     +        mdata,mlat,mlon,nxm,nym,12)
c.diag      filen='airta.dat'                    !Dump to tec.
c.diag      CALL tec1(nxm,nym,12,filen,mlon,mlat,mdata,273.15)
c.diag      filen='airtcoads.dat'                    !Dump to tec.
c.diag      CALL tecij(nxd,nyd,12,filen,dd,273.15)
 
      RETURN
      END
