      SUBROUTINE rNODC(salt,dlon,dlat,zl,nxd,nyd,nzd,nn)
c ------------------------------------------------------------------
c --- Author Knud Simonsen, NERSC 28.3  1996
c ------------------------------------------------------------------
c --- Reads the NODC data and convert them to physical units.
c --- If the values not are witin the allowd range, then they are
c --- set equal tp -999.
c --- Levitus data are stored from 0.5 - 359.5 longitude and 
c --- and from -89.5 to 89.5. In the output the data are stored 
c --- from -179.5 to 179.5.
c -----------------------------------------------------------------
c --- Input:
c --- Output:  salt: salinities    in permill
c ---          zl:   depths in Levitus data
c ---          dlat: Latitudes of the data
c ---          dlon: Longitudes of the data
c ------------------------------------------------------------------
      REAL salt(nxd,nyd,nzd),zl(nzd),
     +     dlat(nyd),dlon(nxd),
     +     smin,smax,tmin,tmax,dums,dumt,dum,
     +     dsaln(10)
      INTEGER i,j,x ,l,nxd2
     

c                                               !Fill in zl
c      DATA zl/0.,10.,20.,30.,50.,75.,100.,125.,150.,200.,250.,
c     +             300.,400.,500.,600.,700.,800.,900.,1000.,1100.,
c     +             1200.,1300.,1400.,1500.,1750.,2000.,2500.,3000.,
c     +             3500.,4000.,4500.,5000.,5500./
      dums= 0. 
      dumt= 10.
      DO i=1,nzd
       IF(i.GE.4)  dumt= 20.  
       IF(i.GE.5)  dumt= 25.
       IF(i.GE.9)  dumt= 50.
       IF(i.GE.12) dumt= 100.
       IF(i.GE.24) dumt = 250.
       IF(i.GE.26) dumt=500.
       zl(i) = dums  
       dums = dums + dumt
      ENDDO

      DO j=1,nyd               !Initialize fields
       DO i=1,nxd
        DO k=1,nzd
         salt(i,j,k) = -999.
        ENDDO
       ENDDO 
      ENDDO
      nxd2 = nxd/2
      smin= 0.         !Minimum allowed salinity
      smax= 100.       !Maximum allowed salinity

      if (nn.eq.1) then
      OPEN(89,FILE='Data/saln_nodc_srf_1x1deg.dat',form='formatted')
      write(*,*) 'reading saln_nodc_srf_1x1deg.dat'
      else
      close(89)
      OPEN(89,FILE='Data/temp_nodc_srf_1x1deg.dat',form='formatted')
      write(*,*) 'reading temp_nodc_srf_1x1deg.dat'
      endif
c
      read(89,'(4i5)') nxf,nyf,nzf,i17
      if ((nxf.ne.nxd).or.(nyf.ne.nyd).or.(nzf.ne.nzd)) then
          write(*,*) nxf,nyf,nzf,nxd,nyd,nzd
          stop '(..)'
      endif
      do k=1,nzf
      read(89,'(10f8.4)') ((salt(i,j,k),i=1,nxf),j=1,nyf)
c
      do i=1,nxd
      do j=1,nyd
      if (salt(i,j,k).lt.-5.) salt(i,j,k)=-999.9999
      enddo
      enddo
c
      enddo
c
      DO i=1,360                 !Longitude 
       dlon(i)=-179.5 + FLOAT(i-1)
      ENDDO
      DO j=1,180                 !Latitude
       dlat(j) = -89.5 + FLOAT(j-1)
      ENDDO

      write(*,*)'NODC data are loaded' 
                        
!     OPEN(3,FILE='salt.dat')
!     write(3,*)'TITLE = "Levitus"'
!     write(3,*)'VARIABLES = "X","Y","S"'
!     write(3,*)'ZONE I= 360, J=180, F=POINT'
!     DO j=1,nyd
!       DO i=1,nxd
!         write(3,'(4f10.2)')dlon(i),dlat(j),salt(i,j,1)
!       ENDDO
!     ENDDO  
!      CLOSE(3)
      GOTO 420 
  400 WRITE(*,*)'Problems to open SALTEMP'
      GOTO 410
  401 WRITE(*,*)'Problems to read data'
  410 STOP 'Sorry'

  420 RETURN
      END
