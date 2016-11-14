      SUBROUTINE rstat(nx,ny,dd,mlon,mlat)
c --- ----------------------------------------------
c --- Purpose: Read data at selected stations from
c --- the file tmpmean.dat and put them into the 
c --- the matrix dd.
c --- ----------------------------------------------
      INTEGER nx,ny,i,j,k,m,stno,id,jd,ifl
      REAL dd(nx,ny,12),ri,rj,mlon(nx,ny),mlat(nx,ny),
     +     lat,lon,stdd(12),fl,dr,drmin
      
      CHARACTER*25 name
      CHARACTER*1  tull
      tull='O'

      OPEN(13,FILE='Data/tmpmean.dat',STATUS='OLD')
      READ(13,*)stno
      WRITE(*,*)'-----------------------------------------------------'
      WRITE(*,'(a6,a18,2a7,a12)')'St. No','     St Name      '
     &       ,'  lat ','  lon  ',' Dist. error'
      WRITE(*,*)'-----------------------------------------------------'

      DO k=1,stno
       READ(13,3,END=99)name,lat, lon
       READ(13,6,END=99)(stdd(m),m=1,12)

       id=1
       jd=1
       drmin=555600.        !App 5 degrees in meters
       DO i=1,nx      
        DO j=1,ny
         dr=spherdist(lon,lat,mlon(i,j),mlat(i,j)) 
         ifl=INT(.5+SIGN(.5,dr-drmin))
         drmin=MIN(dr,drmin) 
         id=ifl*id + (1-ifl)*i 
         jd=ifl*jd + (1-ifl)*j 
        ENDDO
       ENDDO
       IF(drmin.LT.222400.)THEN                 !If within 2 dgr
        DO m=1,12
         
         fl=.5 + SIGN(.5,(dd(id,jd,m)+998.)) ! 1 If contains data, elso 0
         dd(id,jd,m) = fl*dd(id,jd,m) + (1.-fl)*(stdd(m)+273.15)
        ENDDO
        write(53,8)id,jd,name
       ELSE
        write(*,*)'St. at', lon,lat,' is not in model domain'
       ENDIF
       write(*,'(i3,1x,a20,3f7.2,a3)')k,name,lat,lon,drmin*0.001,' km'
      ENDDO
      WRITE(*,*)'----------------------------------------------------'
   99 CLOSE(13) 

    3 FORMAT(A25,2F7.2)
    4 FORMAT(A25,2i4,F7.2)
    6 FORMAT(12f8.2)
    8 FORMAT('TEXT X=',I4,', Y=',I3,
     +    '.,H=0.9,M=GRID,C=YELLOW,T="',A22,' ",BX=FILLED,BXM=0.')

  999 FORMAT(20I4)
  100 FORMAT(10(1x,e10.3))
  101 FORMAT('TITLE="',A10,'"')
  102 FORMAT('VARIABLES="X","Y","lon","lat","CC"')
  103 FORMAT('ZONE T="',A7,'", I=',I3,', J=',I3,', F=BLOCK')
 
      CLOSE(10)
      
      RETURN
      END
