      SUBROUTINE mercator(lon0,glat,glon,mlat,mlon,dx,dy,nx,ny)
c ---------------------------------------------------------------------
c --- Author: Knud Simonsen, NERSC, 01.04.1996 
c --- Modified:
c ---------------------------------------------------------------------
c --- Calc. the distances in rotated Mercator projection of the
c --- positions given by glat,glon relative to the meridiane given
c --- by lon0.
c ---------------------------------------------------------------------
c --- Input: glat,glon: geo. position in scalar point (dgr)
c ---        nx,ny:     dimension of the model grid.
c ---        lon0:      Median,  which the rot. Mercator is
c ---                   realated to. (dgr)
c --- Output:mlat,mlon: Mercator distances (m)
c ---        dx,dy:     (0.1) vector
c ----------------------------------------------------------------------
      INTEGER nx,ny,i,j,mo,im,jm,ip,jp
      REAL  glat(nx,ny),glon(nx,ny) 
     +     ,mlat(nx,ny),mlon(nx,ny)
     +     ,dx(nx,ny),dy(nx,ny)
     +     ,pi,pi2,radian,radinv,fli,flj
     +     ,xlatn,lat,lon,lon0,fl,dum,latd,lond

      data radian/57.29578/,pi/3.14159265/
      pi2 = pi/2.
      radinv=1./radian

c --- ------------------------------------------------------------------------   
c
c.test      lon0=0.                    !Median,  which the rot. Mercator
c.test                                   !Projection is related to.
      DO i=1,nx
       DO j=1,ny
        lon=glon(i,j)-lon0
        fl=.5-SIGN(.5,lon+180.)    !=1 if lon<-180
        lon=lon+fl*360.            !  add 360 
        fl=.5+SIGN(.5,lon-180.)    !=1 if lon>180 
        lon=lon-fl*360.            !  subtract 360
        lat=glat(i,j)

        fl=SIGN(1.,ABS(lon)-90)    !=1  if ABS(lon)>90
                                   !=-1 if ABS(lon)<90
        lon=lon*radinv             !Convert to radians
        lat=lat*radinv

        dum=-SIN(lon)*COS(lat)
        dum=MAX(-1.,MIN(1.,dum)) 
        xlatn=ASIN(dum)              
                                   !Lond dist (in rad) 
                                   !in Mercator projection  
        lond=-ALOG(TAN((2.*xlatn+pi)*.25))

        dum=SIN(lat)/COS(xlatn)
     
        dum=MAX(-1.,MIN(1.,dum))
        latd=ASIN(dum)             !Long dist (in rad) 
        latd=.5*(fl+1.)*pi-fl*latd !in Mercator projection 

        mlon(i,j)=lond*radian*111.2
        mlat(i,j)=latd*radian*111.2


        dy(i,j)=SIN(lon)*SIN(latd)  !(0,1) vector in geo. space
        dx(i,j)=COS(lon)/COS(xlatn) !expressed in Mercator projection
        
       ENDDO
      ENDDO
   99 FORMAT(20I4)
  100 FORMAT(10(1x,e12.6))


c --- ------------------------------------------------------------
c
c      m=1                           !Exp. how to rotate a vector
c      DO i=1,nx
c       DO j=1,ny
c                                    !ud,vd is vector in geospace
c                                    !Angle relative to (1,0) in geospace
c        v1=ATAN2(vd(i,j,m),ud(i,j,m))   
c                                    !Length of vector
c        dl=SQRT(vd(i,j,m)*vd(i,j,m)+ud(i,j,m)*ud(i,j,m))
c                                    !Angle of (1,0) geospace vector 
c                                    !in Mercator projection
c        v2=ATAN2(dy(i,j),dx(i,j))
c                                    !Angle in Mercator space
c        v2=v2+v1
c        ud(i,j,m)=dl*COS(v2) !u in Mercator
c       vd(i,j,m)=dl*SIN(v2) !v in Mercator
c       ENDDO
c      ENDDO
c      write(*,*)'ud,vd rotated from geospace to Mercator'
c      WRITE(1,100)((ud(i,j,1),i=1,nx),j=1,ny) 
c      WRITE(1,100)((vd(i,j,1),i=1,nx),j=1,ny) 
c

cc      stop
c --- ------------------------------------------------------------
      RETURN
      END    


          

