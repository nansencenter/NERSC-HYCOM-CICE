      SUBROUTINE recmwv(du,dv,dlat,dlon,nxd,nyd,
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
c --- Ouput: du,dv: data array in physical units.
c ---        dlon: contains the longitudes.
c ---        dlat: contains the lattitudes
c ----------------------------------------------------------------
      REAL du(nxd,nyd,12),dv(nxd,nyd,12),
     +     dlat(nyd),dlon(nxd)
      CHARACTER filen*15
      REAL a,b,vmin,vmax,valu,valv,
     +     drrho,tabs,tabsmin,tu,tv 
      INTEGER x,y,k,mo ,iu(144),iv(144)

      drrho=1025.*1.5E-6            !rho*dragcoeff.
      tabsmin=drrho*drrho*0.01      !Min Tabs (.1m/sec)
      IF(nxd.NE.144)STOP 'Check  recmwv'

      WRITE(*,*)'Read datafile: ',filen
      OPEN(89,FILE=filen,ACCESS='SEQUENTIAL',FORM='FORMATTED') 

      DO  mo = 1,12                 !For all 12 months do:
       DO  y=1,nyd                  !For all lat. do:
                                    ! (y=1 == 89S
        READ(89,600)(iu(x),x=1,nxd) !Read all values at all long.
        READ(89,600)(iv(x),x=1,nxd) !(x=1 == GRW, read eastw.)
        DO  x=1,nxd                 !Store from 1E to 359E 
                                    !Read wind stresses


         tu = FLOAT(iu(x))/a + b    !Calc.tu
         tv = FLOAT(iv(x))/a + b    !Calc tv.
         tabs=drrho*SQRT(tu*tu+tv*tv)
         tabs=MAX(tabs,tabsmin)
         tabs = SQRT(tabs)
         valu = tu/tabs             !Wind velocity
         valv = tv/tabs            
                                    !If not in accepted range. 
         IF(valu.LT.vmin.OR.valu.GT.vmax)valu = -999. 
         IF(valv.LT.vmin.OR.valv.GT.vmax)valu = -999. 
         du(x,y,mo)= valu 
         dv(x,y,mo)= valv 
        ENDDO  ! x-loop

                                    !Now the data are rearanged.
                                    !Now we start from -179E and
        DO x=1,72                   !end at 179E. 
         k=x+72                     
         valu=  du(x,y,mo)           !Shift u-values
         du(x,y,mo) = du(k,y,mo)
         du(k,y,mo) = valu
         valv=  dv(x,y,mo)           !Shift v-values
         dv(x,y,mo) = dv(k,y,mo)
         dv(k,y,mo) = valv
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
