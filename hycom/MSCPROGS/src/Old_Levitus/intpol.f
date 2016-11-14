      SUBROUTINE intpol(dd,dlat,dlon,nxd,nyd,
     +        mdata,mlat,mlon,nxm,nym,nt)
c ----------------------------------------------------------
c --- This routine interpolates the data from any grid onto
c --- a cartesian grid.
c ----------------------------------------------------------
c --- Author: Knud Simonsen, NERSC, 28.03.1996
c ---
c --- Input: dd: Data file
c ---        dlat: Latitude of the data gridcells.
c ---        dlon: Longitude of the data gridcells. 
c ---        mlat: Latitudes of the model grid.
c ---        mlon: Longitudes of the model gridcells
c ---        nxm,nym: Size of the model grid.
c ---        nxd,nyd: Size of the data grid.
c ---        nt   : number of fields (fi 12 months =12)
c ---
c --- Output: mdata: Data on model grid. 
c ----------------------------------------------------------
c ---
c --- Weighing of data is dw=1/dwa^x, where dwa=MAX(r/drw,distmin)
c --- and x=2 currently. If (dwr/r)^x<dwsmin then no data is included.
c ---
c --- Here drw=1 degree and distmin=1%, i.e. maximum weight
c --- is obtained if data is within 1/100 of a degrees of 
c --- the data point. dwsmin=10% and x=2 means that the influence are
c --- is 1/sqrt(dwsmin) = 3.1 times drw 
c --- 

      PARAMETER(drw=111200.,distmin=0.01,dwsmin=0.1)

c --- ------------------------------------------------------
      REAL mdata(nxm,nym,nt),mlat(nxm,nym),mlon(nxm,nym),
     +          dd(nxd,nyd,nt),dlat(nyd),dlon(nxd)

      INTEGER id,jd,i,j,nxm,nym,nxd,nyd,nt,mo
      REAL yds,ydn,xdw,xde,dr(4),dp(4),dw(4),
     +     dws,dps,ymsm,a1,drwa,dfl
      real, external :: spherdist
      
                     !yd=1 is usually in Antarctic and therefor
                     !we starts from yd=2. If not then put warning.
      IF (dlat(1).GT.-80.) THEN
        WRITE(*,*)'Warning: Check data around Antarctica'
      ENDIF

      mo = 1                          !Fill model file with 
      j = 1                           !no data values which
      DO i=1,nxm*nym*nt               !are equal to  -999.
         mdata(i,j,mo) = -999.  
      ENDDO
                                      !Find the souther most lat.
      ymsm = -85.                      !in the model grid:
!      DO i=1,nxm
!        IF(mlat(i,1).LT.ymsm)THEN
!         id=i
!         ymsm = mlat(i,1)
!        ENDIF
!      ENDDO 
!      ymsm=ymsm-(mlat(id,2)-ymsm+2.)  !Southern limit for search
!      print *,'ymsm=', ymsm
!      stop

                                        
      DO 4 jd=2,nyd                   !For all j in data grid do
       yds = dlat(jd-1)               !Data point south
       ydn = dlat(jd)                 !Data point north
       IF (ymsm.GT.ydn) GOTO 4        !If south of the MD,
                                      !then jump.
       xdw = dlon(nxd)                !First data point west
       DO 5 id=1,nxd                  !For all i in data grid do  
        xde = dlon(id)                !Data point east
        DO 6 j=1,nym                  !For all (i,j) in the model 
         DO 7 i=1,nxm                 !grid find the model grid
                                      !cell, which is between the
                                      !data points. 

                                      !Check model latitude
          IF (mlat(i,j).GT.ydn.OR.mlat(i,j).LE.yds) GOTO 7
           
          IF (id.EQ.1) THEN           !In first point ~180W special
                                      !tratment is needed.
            IF (ABS(mlon(i,j)).LE.xdw.OR.
     +          ABS(mlon(i,j)).LT.xde) GOTO 7
          ELSE 
            IF (mlon(i,j).LE.xdw.OR.mlon(i,j).GT.xde) GOTO 7
          ENDIF 

          dr(1)=spherdist(xdw,ydn,mlon(i,j),mlat(i,j))
          dr(2)=spherdist(xdw,yds,mlon(i,j),mlat(i,j))
          dr(3)=spherdist(xde,ydn,mlon(i,j),mlat(i,j))
          dr(4)=spherdist(xde,yds,mlon(i,j),mlat(i,j))

c.diag          write(*,*)id,jd,'lat',mlat(i,j),ydn,yds
c.diag          write(*,*)i,j,'lon',mlon(i,j),xdw,xde
c.diag          write(*,*)dr(1),dr(2),dr(3),dr(4)


          DO mo=1,nt                        !For all months do  
            dp(1) = dd(id-1,jd,mo)          !NW datavalue
            dp(2) = dd(id-1,jd-1,mo)        !SW datavalue
            dp(3) = dd(id,jd,mo)            !NE datavalue
            dp(4) = dd(id,jd-1,mo)          !SE datavalue
            dws = 0.
            dps = 0.
            DO k=1,4 
                                            !If no data (=-999) then
                                            !weight equals zero, otherwise
                                            !inverse distance
              drwa=MAX(dr(k)/drw,distmin)
              drwa=drwa*drwa

              dw(k)=(.5+SIGN(.5,998.+dp(k)))/drwa
              dws = dws + dw(k)             !Sum up weights
              dps = dps + dw(k)*dp(k)       !Sum up weighted data
            ENDDO  
            dfl=.5+SIGN(.5,dws-dwsmin)      !=1 if dws=>dwsmin, =>mdata=dps/dws
                                            !=0 if dws<dwsmin,  =>mdata=-999 
            mdata(i,j,mo)=dfl*dps/(dws+1.-dfl) - (1.-dfl)*999.
c.diag            write(*,*)dfl,dps,dws,mdata(i,j,mo)
          ENDDO
           
   7     CONTINUE 
   6    CONTINUE
        xdw=xde                       !Shift data point east to west
   5   CONTINUE 
   4  CONTINUE
      RETURN
      END 
