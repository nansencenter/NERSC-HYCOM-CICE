      SUBROUTINE rrain(dd,dlat,dlon,nxd,nyd,
     +               filen,a,b,vmin,vmax)
c -----------------------------------------------------------------
c --- Author Knud Simonsen, NERSC,2/12 1993.
c --- -------------------------------------------------------------
c --- Reads the Legates and Willmott precipitation data
c ---  data and converts to physical units (SI units)
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
      INTEGER x,y,k,mo ,id(360)

      IF(nxd.NE.360)STOP 'Check rrain'

      WRITE(*,*)'Read datafile: ',filen
      OPEN(89,FILE=filen,ACCESS='SEQUENTIAL',FORM='FORMATTED',ERR=602) 

      DO  mo = 1,12                 !For all 12 months do:
       DO  y=1,nyd                  !For all lat. do:
                                    ! (y=1 == 89S
        READ(89,600)(id(x),x=1,360) !Read all values at all long.
                                    !(x=1 == GRW, read eastw.)
        DO  x=1,360                 !Store from 1E to 359E 
         val = FLOAT(id(x))/a + b   !Calc. valu
         IF(val.LT.vmin.OR.val.GT.vmax) THEN
          val = -999.               !If not in accepted range.
         ENDIF
         dd(x,y,mo)= val 
        ENDDO  ! x-loop
                                    !Now the data are rearanged.
                                    !Now we start from -179E and
        DO x=1,180                   !end at 179E. 
         k=x+180                     
         val=  dd(x,y,mo)           !Shift u-values
         dd(x,y,mo) = dd(k,y,mo)
         dd(k,y,mo) = val
        ENDDO   !x-loop
       ENDDO   !y-loop
      ENDDO   ! mo-loop   
      CLOSE(89)
      DO y=1,nyd                     !Calc. lattitude
        dlat(y)= -89.5+ FLOAT(y-1)
      ENDDO
      DO x=1,360                     !Calc. longitude.
        dlon(x) = -179.5 + FLOAT(x-1)
      ENDDO
  600 FORMAT(20I4)
      GOTO 603
  602 STOP 'File Error'
  603 RETURN
      END
