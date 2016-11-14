      SUBROUTINE recmwf(dd,dlat,dlon,nxd,nyd,
     +               filen,a,b,vmin,vmax)
c -----------------------------------------------------------------
c --- Author Knud Simonsen, NERSC,12/01 1994.
c --- -------------------------------------------------------------
c --- Reads the ECMWF data and converts to physical units (SI units)
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
      INTEGER x,y,k,mo ,id(144)

      IF(nxd.NE.144)STOP 'Check rcoads'

      WRITE(*,*)'Read datafile: ',filen
      OPEN(89,FILE=filen,ACCESS='SEQUENTIAL',FORM='FORMATTED') 

      DO  mo = 1,12                 !For all 12 months do:
       DO  y=1,nyd                  !For all lat. do:
                                    ! (y=1 == 89S
        READ(89,600)(id(x),x=1,nxd) !Read all values at all long.
                                    !(x=1 == GRW, read eastw.)
        DO  x=1,nxd                 !Store from 1E to 359E 
         val = FLOAT(id(x))/a + b   !Calc. valu
         IF(val.LT.vmin.OR.val.GT.vmax.OR.y.EQ.nyd) THEN
          val = -999.               !If not in accepted range.
         ENDIF
         dd(x,y,mo)= val 
        ENDDO  ! x-loop
                                    !Now the data are rearanged.
                                    !Now we start from -179E and
        DO x=1,72                   !end at 179E. 
         k=x+72                     
         val=  dd(x,y,mo)           !Shift u-values
         dd(x,y,mo) = dd(k,y,mo)
         dd(k,y,mo) = val
        ENDDO   !x-loop
       ENDDO   !y-loop
      ENDDO   ! mo-loop   
      CLOSE(89)
      DO y=1,nyd                     !Calc. lattitude
        dlat(y)= FLOAT(y-37)*2.5
      ENDDO
      DO x=1,nxd                     !Calc. longitude.
        dlon(x) = FLOAT(x-72)*2.5
      ENDDO
  600 FORMAT(20I4)
      RETURN
      END
