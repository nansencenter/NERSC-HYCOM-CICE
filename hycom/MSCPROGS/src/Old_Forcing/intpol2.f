      SUBROUTINE intpol2(dd,dlat,dlon,nxd,nyd,
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
c --- and x=2 currently. If (dwr/r)<dwsmin then no data is included.
c ---
c --- Here drw=1 degree and distmin=1%, i.e. maximum weight
c --- is obtained if data is within 1/100 of a degrees of 
c --- the data point. dwsmin=10% and x=2 means that the influence are
c --- is 1/sqrt(dwsmin) = 3.1 times drw 
c --- 

      PARAMETER(drw=111200.,distmin=0.01,dwsmin=0.1)
      PARAMETER(nstoc=4)

c --- ------------------------------------------------------
      INTEGER id,jd,i,j,nxm,nym,nxd,nyd,nt,mo
     +       ,idd(nstoc),jdd(nstoc)

      REAL mdata(nxm,nym,nt),mlat(nxm,nym),mlon(nxm,nym),
     +          dd(nxd,nyd,nt),dlat(nxd,nyd),dlon(nxd,nyd)


      REAL dr(nstoc),drr,dw(nstoc),dp(nstoc)

      LOGICAL dflg


 
      mo = 1                          !Fill model file with 
      j = 1                           !no data values which
      DO i=1,nxm*nym*nt               !are equal to  -999.
         mdata(i,j,mo) = -999.  
      ENDDO

      DO  j=1,nym                  !For all (i,j) in the model
      write(*,*) j
       DO  i=1,nxm                 !grid find the model grid

c      DO  j=50,50                  !For all (i,j) in the model
c       DO  i=110,110                 !grid find the model grid
        
        DO k=1,nstoc
         dr(k)=1.E29
         idd(k)=0
         jdd(k)=0
        ENDDO  
                                         
        DO jd=1,nyd                !For all j in data grid do
         DO id=1,nxd               !For all i in data grid do  
                                 
                                   !Compute disstance between
                                   !modelpint and datapoints
     
                                    
          drr=spherdist(dlon(id  ,jd  ),dlat(id  ,jd  )
     &                    ,mlon(i,j),mlat(i,j))
 
          dflg=.true. 
          DO k=1,nstoc,1

           IF(drr.LT.dr(k).AND.dflg)THEN
            IF(k.LT.nstoc)THEN      !Restoc
             DO k1=nstoc,k+1,-1
              dr(k1)=dr(k1-1)
              idd(k1)=idd(k1-1)
              jdd(k1)=jdd(k1-1)
             ENDDO
            ENDIF
            dr(k)=drr
            idd(k)=id
            jdd(k)=jd
            dflg=.false.
           ENDIF 
          ENDDO
         ENDDO  !id
        ENDDO !jd
c        write(*,*)dr,idd,jdd



        DO mo=1,nt                        !For all months do  
         dfl=0.
         dws=0.
         dps=0.
         DO k=1,nstoc
          dp(k)=dd(idd(k),jdd(k),mo)
                                          !If no data (=-999) then
                                          !weight equals zero, otherwise
                                          !inverse distance
          drwa=MAX(dr(k),distmin)
          drwa=1./(drwa*drwa)
          
c          dfl=dfl+(.5+SIGN(.5,dp(k)+0.001))  !Number of OK-data
          dfl=dlf+(.5+SIGN(.5,dp(k)+0.001))   !Number of OK data          
          dw(k)=(.5+SIGN(.5,dp(k)+0.001))*drwa
          dws = dws + dw(k)             !Sum up weights
          dps = dps + dw(k)*dp(k)       !Sum up weighted data
         ENDDO   !k
      
         dfl=.5+SIGN(.5,dfl-0.5)         !=1 if OK 
            
         mdata(i,j,mo)=dfl*dps/(dws+1.-dfl) - (1.-dfl)*999.
        ENDDO !mo
       ENDDO
      ENDDO

      RETURN
      END 
