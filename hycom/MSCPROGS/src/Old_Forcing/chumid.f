      SUBROUTINE chumid(humid,iw,nx,ny,mo,
     +        iflg,mlon,mlat)

      INTEGER iw,nx,ny,mo,iflg(nx,ny)
      REAL humid(nx,ny,mo),mlon(nx,ny),mlat(nx,ny),
     +    a,b,vmin,vmax,ncc(12)
      CHARACTER*15 filen
       
     
      IF(iw.EQ.1.OR.iw.EQ.2) THEN 
       filen= 'Data/WETNESS'         !Read Relative Humidity
       a= 1000.
       b= 0.
       vmin=0.1
       vmax=1.0
       CALL lcoads(humid,iflg,nx,ny,mlon,mlat,filen,a,b,
     +            vmin,vmax)
       IF(iw.EQ.2) THEN                 !Include values from Maykut, 1978,
                                        !JGR, Vol. 83, 3646-3658.,Table 1.
        ncc(1) = .94       
        ncc(2) = .91      
        ncc(3) = .90
        ncc(4) = .90
        ncc(5) = .90
        ncc(6) = .90
        ncc(7) = .90
        ncc(8) = .90
        ncc(9) = .90
        ncc(10)= .90
        ncc(11)= .91
        ncc(12)= .96

        write(*,*)'-----------------------------'
        write(*,*)'Includes RH from Maykut, 1978'
        write(*,*)'See subroutine chumid.f      '
        write(*,*)'-----------------------------'
     

        DO i=1,nx
         DO j=1,ny
          IF(mlat(i,j).GE.85.)THEN
           IF(mlon(i,j).LT.-90..OR.mlon(i,j).GT.90.)THEN
            DO k=1,12
             fl = .5+SIGN(.5,(998.+humid(i,j,k)))     != 0 if no data
c.diag             write(*,*)humid(i,j,k),fl,ncc(k)
             humid(i,j,k)=humid(i,j,k)*fl + (1.-fl)*ncc(k)
            ENDDO
           ENDIF
          ENDIF
         ENDDO
        ENDDO
       ENDIF !iw=2
                                      !Fillup missing data
                                      !and land mask.
       CALL fillup(humid,nx,ny,iflg,12,0.8)

      ELSE
       WRITE(*,*)'Rel.Humid. does not exist'
      ENDIF
      RETURN
      END



