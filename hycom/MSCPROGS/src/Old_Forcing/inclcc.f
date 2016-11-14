      SUBROUTINE inclcc(cloud,nx,ny,mlon,mlat)
c --- ----------------------------------------------------
c --- Insert cloud cover data from Husche, 1969 as cited by
c --- Holland et al, 1993 in the cloud cover filed in the 
c --- area north of 85N between 90 -- -90E (90, 180, -180,-90)
c --- ----------------------------------------------------
 
      INTEGER nx,ny
      REAL cloud(nx,ny,12) ,
     +   mlon(nx,ny),mlat(nx,ny),
     +   ncc(12),ri,rj,lat,lon,fl


      ncc(1) = .50                         !Cloud cover in the central Arctic
      ncc(2) = .50                         !from  Holland et al, 1993.
      ncc(3) = .50
      ncc(4) = .55
      ncc(5) = .70
      ncc(6) = .75
      ncc(7) = .75
      ncc(8) = .80
      ncc(9) = .80
      ncc(10)= .70
      ncc(11)= .60
      ncc(12)= .50
      write(*,*)'----------------------------'
      write(*,*)'Include CC from Husche, 1969'
      write(*,*)'See subroutine inclcc.f     '
      write(*,*)'----------------------------'
      DO i=1,nx
       DO j=1,ny
        IF(mlat(i,j).GE.85.)THEN
         IF(mlon(i,j).LT.-90..OR.mlon(i,j).GT.90.)THEN
          DO k=1,12
           fl = .5+SIGN(.5,(998.+cloud(i,j,k)))
c.diag           write(*,*)cloud(i,j,k),fl,ncc(k)
           cloud(i,j,k)=cloud(i,j,k)*fl + (1.-fl)*ncc(k)
          ENDDO 
         ENDIF
        ENDIF 
       ENDDO
      ENDDO
      RETURN
      END



