      SUBROUTINE rcoads(dd,dlat,dlon,nxd,nyd,
     +               filen,a,b,vmin,vmax)
c -----------------------------------------------------------------
c --- Author Knud Simonsen, NERSC,2/12 1993.
c --- -------------------------------------------------------------
c --- Reads the COADS data and converts to physical units (SI units)
c -----------------------------------------------------------------
c --- Input: filen: filaname 
c ---        nxd,nyd: size of data file. 
c ---        a,b: const. to convert.
c ---        vmin,vmax: min and max accepted values
c --- Ouput: dd: data array in physical units.
c ---        dlon: contains the longitudes.
c ---        dlat: contains the lattitudes
c ----------------------------------------------------------------
      REAL dd(nxd,nyd,12),
     +     dlat(nyd),dlon(nxd)
      CHARACTER filen*15
      REAL a,b,vmin,vmax,val,td
      INTEGER x,y,k,mo ,id(180)

      IF(nxd.NE.180)STOP 'Check rcoads'

      WRITE(*,*)'Read datafile: ',filen
      OPEN(89,FILE=filen,ACCESS='SEQUENTIAL',FORM='FORMATTED',ERR=602) 

      DO  mo = 1,12                 !For all 12 months do:
       DO  y=1,nyd                  !For all lat. do:
                                    ! (y=1 == 89S
        READ(89,600)(id(x),x=1,180) !Read all values at all long.
                                    !(x=1 == GRW, read eastw.)
        DO  x=1,180                 !Store from 1E to 359E 
         val = FLOAT(id(x))/a + b   !Calc. valu
         IF(val.LT.vmin.OR.val.GT.vmax) THEN
          val = -999.               !If not in accepted range.
         ENDIF
         dd(x,y,mo)= val 
        ENDDO  ! x-loop
                                    !Now the data are rearanged.
                                    !Now we start from -179E and
        DO x=1,90                   !end at 179E. 
         k=x+90                     
         val=  dd(x,y,mo)           !Shift u-values
         dd(x,y,mo) = dd(k,y,mo)
         dd(k,y,mo) = val
        ENDDO   !x-loop
       ENDDO   !y-loop
c.ks       write(*,'(8f8.2)')dd(87,70,mo)
c.ks     &           ,dd(97,60,mo),dd(90,58,mo),dd(96,78,mo)
c.ks     &           ,dd(87,75,mo),dd(90,63,mo),dd(88,58,mo)
c.ks     &           ,dd(91,68,mo)
      ENDDO   ! mo-loop   
      CLOSE(89)
      DO y=1,nyd                     !Calc. lattitude
        dlat(y)= -FLOAT(nyd-1)+ FLOAT(y-1)*2.
      ENDDO
      DO x=1,180                     !Calc. longitude.
        dlon(x) = -179. + FLOAT(x-1)*2.
      ENDDO
  600 FORMAT(20I4)
      GOTO 603
  602 STOP 'File Error'
  603 RETURN
      END
