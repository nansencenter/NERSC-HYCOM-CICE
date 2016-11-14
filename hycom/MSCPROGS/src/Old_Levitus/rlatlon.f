      SUBROUTINE rlatlon(mlat,mlon)
c ----------------------------------------------------------------
c --- Author Knud Simonsen, 28-03 1996, NERSC
c ---  partly adopted from SPEM grid generation package, developed 
c --- by Kate Hedstrom and ohn Wilkin 
c ---
c --- Modified by HD 19.09.1996
c ----------------------------------------------------------------
c --- Load geo. pos of the model grid at scalar and u, v points
c ----------------------------------------------------------------
c --- Input: L,M    !Model size
c --- Output:  
c ----------------------------------------------------------------
c.test      PARAMETER(L=101,M=81)

      INCLUDE 'dimension.h'

      INTEGER nxl,nyl
      CHARACTER*80    gridid

      REAL platl(0:idm+1,0:jdm+1),plonl(0:idm+1,0:jdm+1)           !pos in pressure point
     .    ,qlatl(0:idm+1,0:jdm+1),qlonl(0:idm+1,0:jdm+1)           !pos in vorticity point
     .    ,ulatl(0:idm+1,0:jdm+1),ulonl(0:idm+1,0:jdm+1)           !pos in u-points
     .    ,vlatl(0:idm+1,0:jdm+1),vlonl(0:idm+1,0:jdm+1)           !pos in v-points
     .    ,mlat(idm,jdm),mlon(idm,jdm)
      

      WRITE(*,*)'Load grid positions from file: latlon.dat'
c
      open(73,file='Data/latlon.dat',form='formatted')
      READ(73,'(2i5)',ERR=250)nxl,nyl
      write(*,*)nxl,nyl
      IF(nxl.NE.idm.OR.nyl.NE.jdm)THEN
        WRITE(*,*)'Mismatch between lat-lon file dimension'
        WRITE(*,*)'and model dimension. I quit!!!!'
        CLOSE(73)
        STOP'Error'
      ENDIF
c
      read(73,'(a80)')gridid
      read(73,300)((platl(i,j),i=1,idm),j=1,jdm)
      read(73,300)((plonl(i,j),i=1,idm),j=1,jdm)
      read(73,300)((qlatl(i,j),i=0,idm),j=0,jdm)
      read(73,300)((qlonl(i,j),i=0,idm),j=0,jdm)
      read(73,300)((ulatl(i,j),i=1,idm),j=0,jdm)
      read(73,300)((ulonl(i,j),i=1,idm),j=0,jdm)
      read(73,300)((vlatl(i,j),i=0,idm),j=1,jdm)
      read(73,300)((vlonl(i,j),i=0,idm),j=1,jdm)
      close(73)

 300  FORMAT(15e14.7)

      DO i=1,idm
       DO j=1,jdm
        mlon(i,j)=plonl(i,j)
        mlat(i,j)=platl(i,j)
       ENDDO
      ENDDO



c.ks      WRITE(*,*)'TITLE=""'
c.ks      WRITE(*,*)'VARIABLES="psilat","psilon","rholat","rholon",',
c.ks     &           '"ulat","ulon","vlat","vlon"'
c.ks      WRITE(*,*)'ZONE  I=',L,', J=',M,', F=BLOCK'
c.ks      WRITE(*,800)((psilat(i,j), i=1,L),j=1,M)
c.ks      WRITE(*,800)((psilon(i,j), i=1,L),j=1,M)
c.ks      WRITE(*,800)((rholat(i,j), i=1,L),j=1,M)
c.ks      WRITE(*,800)((rholon(i,j), i=1,L),j=1,M)
c.ks      WRITE(*,800)((ulat(i,j),   i=1,L),j=1,M)
c.ks      WRITE(*,800)((ulon(i,j),   i=1,L),j=1,M)
c.ks      WRITE(*,800)((vlat(i,j),   i=1,L),j=1,M)
c.ks      WRITE(*,800)((vlon(i,j),   i=1,L),j=1,M)
 800  FORMAT(20f8.2)

      GOTO 900
 250  WRITE(*,*)'I can not read the file, Sorry'
      CLOSE(73)
      STOP'Error'
 900  RETURN
      END 
    
