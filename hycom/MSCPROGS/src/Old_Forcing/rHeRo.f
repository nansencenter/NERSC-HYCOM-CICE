      SUBROUTINE rHeRo(du,dv,dlat,dlon,nxd,nyd)
c -----------------------------------------------------------------
c --- Author Knud Simonsen, NERSC,2/12 1993.
c --- -------------------------------------------------------------
c --- Reads the psuedo wind stresses from Hellerman and Rosenstein 
c --- and converts to physical units. The stored pseudo stresses
c --- are tx=sqrt(c*rho)*u, ty=sqrt(c*rho)*v, where c=1.5E-6 is
c --- the atmospheric dragcoefficient, rho=1025 a density and (u,v)
c --- the wind velocity.
c --- Ref.:Subroutine WIND in OPYC, p. 89 in OPYC-manual.
c --- The values are checked if they are  within allowed physical 
c --- range. If not in the allowed range, then the value is set equal 
c --- to -999. 
c -----------------------------------------------------------------
c --- Input: 
c ---        nxd,nyd: size of data file. 
c --- Ouput: du: u-component of the windsress.
c ---        dv: v-component of the windsress.
c ---        dlon: contains the longitudes.
c ---        dlat: contains the lattitudes
c ----------------------------------------------------------------
      REAL du(nxd,nyd,12),dv(nxd,nyd,12),
     +     dlat(nyd),dlon(nxd)
      CHARACTER filen*15
      REAL a,b,vmin,vmax,valu,valv,tu,tv,drrho
      INTEGER x,y,k,mo,iu(180),iv(180),nxd,nyd

      IF(nxd.NE.180.OR.nyd.NE.90)STOP 'Check rHeRo'

      filen='Data/UVGLOB'
      a=5000.
      drrho=SQRT(1025.*1.5E-6)      !SQRT(rho*dragcoeff.)
      vmin = -  30.                 !Min. value (-20 m/s) 
      vmax =    30.                 !Max. value (20m/s)

      WRITE(*,*)'Read datafile: ',filen
      OPEN(89,FILE=filen,ACCESS='SEQUENTIAL',FORM='FORMATTED') 
                                  !The first block of data contains 
      READ(89,600)(iu(x),x=1,nxd) !tha annual mean. These data are
      READ(89,600)(iv(x),x=1,nxd) !dumped.

      DO  mo = 1,12                 !For all 12 months do:
       DO  y=1,nyd                  !For all lat. do:
                                    ! (y=1 == 89S
        READ(89,600)(iu(x),x=1,nxd) !Read all values at all long.
        READ(89,600)(iv(x),x=1,nxd) !(x=1 == GRW, read eastw.)
        DO  x=1,nxd                 !Store from 1E to 359E 
         valu = (iu(x)-a) / a       !Calc. u-comp
         valv = (iv(x)-a) / a       !Calc. v-comp
c         tu = valu*SQRT(valu*valu+valv*valv) !Wind stresses
c         tv = valv*SQRT(valu*valu+valv*valv)
         tu =valu/drrho             !Wind velocities
         tv =valv/drrho
         IF(tu.LT.vmin.OR.tu.GT.vmax) THEN
          tu = -999.               !If not in accepted range.
         ENDIF
         IF(tv.LT.vmin.OR.tv.GT.vmax) THEN
          tv = -999.               !If not in accepted range.
         ENDIF
         du(x,y,mo)= tu
         dv(x,y,mo)= tv
        ENDDO  ! x-loop
                                    !Now the data are rearanged.
                                    !Now we start from -179E and
        DO x=1,90                   !end at 179E. 
         k=x+90                     
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
      DO y=1,90                      !Calc. lattitude
        dlat(y)= -89. + FLOAT(y-1)*2.
      ENDDO
      DO x=1,180                     !Calc. longitude.
        dlon(x) = -179. + FLOAT(x-1)*2.
      ENDDO
      write(*,*)'ok rHeRo.f'
  600 FORMAT(20I4)
      RETURN
      END
