       program  protorot
c       PARAMETER (nx=121,ny=71,nz=1)
       PARAMETER (nx=181,ny=240,nz=1)

c ---------------------------------------------------------------------
c --- Author: Knud Simonsen, NERSC, 01.04.1996 
c --- Modified:
c ---------------------------------------------------------------------
c --- Rotates a vector field from a grid in geographical coordinates
c --- into the grid defined by the lat,lon in the input variables.
c --- C-grid assumed.
c --- 
c --- Below is also included rutines which convert from model grid
c --- to geographical coordinates, and from geographical coordinates
c --- into a rotated Mercator projection. These rutines are extracted
c --- from this subroutine into another routine  
c ---------------------------------------------------------------------
c --- Input: mlat, mlon: position in scalar point
c ---        ulat,ulon: position in u-point
c ---        vlat,vlon: position in v-point
c ---        nx,ny:     dimension of the model grid.
c ---        ud,vd:  Unrotated vector components. 
c --- Output: ud,vd: Rotated vector components.
c ----------------------------------------------------------------------
      INTEGER nx,ny,i,j,mo,im,jm,ip,jp,nxl,nyl
      REAL ud(nx,ny,nz),vd(nx,ny,nz)
     +     ,mlat(nx,ny),mlon(nx,ny)
     +     ,ulat(nx,0:ny),ulon(nx,0:ny)
     +     ,vlat(0:nx,ny),vlon(0:nx,ny)
     +     ,qlat(0:nx,0:ny),qlon(0:nx,0:ny)      !pos in vorticity point
     +     ,pi,pi2,radian,radinv,fli,flj
     +     ,u_up,v_up,u_vp,v_vp,theta_up,theta_vp
     +     ,dlon
     +     ,xlatn ,v1,v2,dl,u_pp,v_pp
     +     ,lat,lon,lon0,fl,dum,latd,lond,dx(nx,ny),dy(nx,ny)
     +     ,dr,dra(4)
     +     ,mlat1(nx,ny),mlon1(nx,ny)

      LOGICAL geo_to_cl,cl_to_geo,geo_to_merc,polcorr


      data radian/57.29578/,pi/3.14159265/
      pi2 = pi/2.
      radinv=1./radian

c --- -----------------------------------------------------------
c --- Test to be performed:
c ---
      geo_to_cl=.FALSE.          !From geospace to curvelinear grid
      cl_to_geo=.TRUE.           !From curvelinear to geospace
      polcorr = .False.          !Correction around the Pole. 
      geo_to_merc=.TRUE.        !From geospace to Mercator     
      
c --
c --- ---------------------------------------- Load positions ---
      OPEN(73,file='latlon.dat',form='formatted')
      READ(73,'(2i5)',ERR=250)nxl,nyl
      write(*,*)nx,ny
      IF(nxl.NE.nx.OR.nyl.NE.ny)THEN
        WRITE(*,*)'Mismatch between lat-lon file dimension'
        WRITE(*,*)'and model dimension. I quit!!!!'
        CLOSE(73)
        STOP'Error'
      ENDIF

      READ(73,'(a80)',ERR=250)gridid
      write(*,*)'gridid',gridid
      READ(73,'(15e14.7)',ERR=250)mlat,mlon   !In pressure point
      READ(73,'(15e14.7)',ERR=250)qlat,qlon   !In vorticity point
      READ(73,'(15e14.7)',ERR=250)ulat,ulon   !In u-points
      READ(73,'(15e14.7)',ERR=250)vlat,vlon   !In v-points
c      READ(73,'(15e14.7)')r3
      CLOSE(73)
c      write(*,*) 'read 73 ok'
      GOTO 251
 250  WRITE(*,*)'Error when readinf file: latlon.dat'
      STOP'Error'
 251  CONTINUE


c --- ---------------------------------------- test wind -------
      DO i=1,nx
       DO j=1,ny
        ud(i,j,1)=1. 
        vd(i,j,1)=0. 
       ENDDO
      ENDDO
c
      write(*,*)'i,j'
      WRITE(1,'(20i4)')((i,i=1,nx),j=1,ny)
      WRITE(1,'(20i4)')((j,i=1,nx),j=1,ny)

      write(*,*)'mlon,mlat, unmodified'
      WRITE(1,100)((mlon(i,j),i=1,nx),j=1,ny)
      WRITE(1,100)((mlat(i,j),i=1,nx),j=1,ny)
      write(*,*)'ud,vd test values'
      WRITE(1,100)((ud(i,j,1),i=1,nx),j=1,ny) 
      WRITE(1,100)((vd(i,j,1),i=1,nx),j=1,ny) 

   99 FORMAT(20I4)
  100 FORMAT(10(1x,e12.6))
c.test      GOTO 900

c --- ----------------------------------------------------------
c --- Assumes that all parameters are provided in scalar point 
c --- and interpolates into the U- and V (C-grid) points, and
c --- perform the rotation rquired in curvlinear grid.
c --- ----------------------------------------------------------                             
      IF(geo_to_cl)THEN
      DO j=ny,1,-1
       DO i=nx,1,-1
        im=MAX(i-1,1)
        fli=FLOAT(i-im)
        jm=MAX(j-1,1)
        flj=FLOAT(j-jm)
                                       !Rotation angle in u-point 
        dlon=mlon(i,j)-mlon(im,j)
        IF(dlon.LT.-180.)dlon=360.0+dlon
        IF(dlon.GT.180.)dlon=dlon-360.0
        theta_up = atan2(mlat(i,j)-mlat(im,j), 
     &            dlon*cos(radinv*.5*(mlat(i,j)+mlat(im,j))) )
c
                                       !Rotation angle in v-point 
        dlon=mlon(i,j)-mlon(i,jm)
        IF(dlon.LT.-180.)dlon=360.0+dlon
        IF(dlon.GT.180.)dlon=dlon-360.0
        theta_vp = atan2(mlat(i,j)-mlat(i,jm),
     &            dlon*cos(radinv*.5*(mlat(i,j)+mlat(i,jm))) )
     &            -pi2

        DO mo=1,nz 
         u_up=.5*(ud(i,j,mo)+ud(im,j,mo)) !Unrotated vel. in u-point
         v_up=.5*(vd(i,j,mo)+vd(im,j,mo))
 
         u_vp=.5*(ud(i,j,mo)+ud(i,jm,mo)) !Unrotated vel. in v-point
         v_vp=.5*(vd(i,j,mo)+vd(i,jm,mo))
                                          !Perform rotation:
         ud(i,j,mo)= (u_up*COS(theta_up)+ v_up*SIN(theta_up))*fli
         vd(i,j,mo)=(-u_up*SIN(theta_vp)+ v_up*COS(theta_vp))*flj

c         ud(i,j,mo)= theta_up
c         vd(i,j,mo)=theta_vp

cc --- ---------------------------------------------------- test --------
c.test         vd(i,j,mo)=(u_vp*SIN(theta_vp)+ v_vp*COS(theta_vp))*flj

c.diag         write(*,*)i,j,theta_up*radian,theta_vp*radian,
c.diag     &      ud(i,j,mo),vd(i,j,mo)
c --- ---------------------------------------------------- test --------
  
        ENDDO  !i-loop
       ENDDO  !j-loop
      END DO  !mo-loop

                                          !Perform interpolation over the
                                          !polar area in order to hide the
                                          !error arising from the calc. of
                                          !the rot. angle

c      im=103
c      ip=109
c      jm=42
c      jp=55 
c      
c      DO i=im+1,ip-1
c       DO j=jm+1,jp-1
c        dra(1)=spherdist(mlon(i,jm),mlat(i,jm),mlon(i,j),mlat(i,j)) 
c        dra(2)=spherdist(mlon(i,jp),mlat(i,jp),mlon(i,j),mlat(i,j)) 
c        dra(3)=spherdist(mlon(im,j),mlat(im,j),mlon(i,j),mlat(i,j)) 
c        dra(4)=spherdist(mlon(ip,j),mlat(ip,j),mlon(i,j),mlat(i,j)) 
c        dum=dra(1)+dra(2)+dra(3)+dra(4)
c        ud(i,j,mo)=(ud(i,jm,mo)*dra(2)+ud(i,jp,mo)*dra(1)
c     &             + ud(im,j,mo)*dra(4)+ud(ip,j,mo)*dra(3))/dum
c        ud(i,j,mo)=(vd(i,jm,mo)*dra(2)+vd(i,jm,mo)*dra(1)
c     &             + vd(im,j,mo)*dra(4)+vd(ip,j,mo)*dra(3))/dum
c         write(*,*)dra,dum
c       ENDDO
c      ENDDO  

 900  CONTINUE
      write(*,*)'ud,vd rotated from geo space to xi,eta space'
      WRITE(1,100)((ud(i,j,1),i=1,nx),j=1,ny)
      WRITE(1,100)((vd(i,j,1),i=1,nx),j=1,ny)
      ENDIF 

c --- ------------------------------------------------------------------------   
c --- Assumes that ud,vd contains velocities in xi, eta space on C-grid
c --- and rotates to geospace. 
c --- ------------------------------------------------------------------------   
       IF(cl_to_geo)THEN
        CALL curlin_geo(ud,vd,mlon,mlat,ulon,ulat,vlon,vlat,nx,ny,nz)

       write(*,*)'ud,vd rotated from xi, eta to   geo. space'
       WRITE(1,100)((ud(i,j,1),i=1,nx),j=1,ny) 
       WRITE(1,100)((vd(i,j,1),i=1,nx),j=1,ny) 
       WRITE(1,100)((ud(i,j,1),i=1,nx),j=1,ny) 
       WRITE(1,100)((vd(i,j,1),i=1,nx),j=1,ny) 
      ENDIF ! cl_to_geo

c --- --------------------------------------------------------------
c --- Assums that  ud,vd are in geospace i.e ud=E direction and
c --- vd=N diection and rotates into Mercator projection
c --- --------------------------------------------------------------
      IF(geo_to_merc)THEN

      lon0=-40.                    !Median,  which the rot. Mercator
                                   !Projection is related to.
      CALL mercator(-40.,mlat,mlon,mlat1,mlon1,dx,dy,nx,ny)    

c      WRITE(1,100)((dx(i,j),i=1,nx),j=1,ny)
c      WRITE(1,100)((dy(i,j),i=1,nx),j=1,ny)

      write(*,*)'mlon,mlat in Mercator Projection'
c      WRITE(1,100)((mlon(i,j),i=1,nx),j=1,ny)
c      WRITE(1,100)((mlat(i,j),i=1,nx),j=1,ny)
      WRITE(1,100)((mlon1(i,j),i=1,nx),j=1,ny)
      WRITE(1,100)((mlat1(i,j),i=1,nx),j=1,ny)
cc      WRITE(1,100)(((mlon(MIN(i+1,i),j)-mlon(i,j)),
cc     &          i=1,nx),j=1,ny)
cc      WRITE(1,100)(((-mlat(i,MIN(j+1,j))+mlat(i,j)),i=1,nx),j=1,ny)
cc      write(*,*)'dx,dy in Mercator Projection'
c      WRITE(1,100)((dx(i,j),i=1,nx),j=1,ny)
c      WRITE(1,100)((dy(i,j),i=1,nx),j=1,ny)
c
      m=1
      DO i=1,nx
       DO j=1,ny
                                    !ud,vd is vector in geospace
                                    !Angle relative to (1,0) in geospace
        v1=ATAN2(vd(i,j,m),ud(i,j,m))   
                                    !Length of vector
        dl=SQRT(vd(i,j,m)*vd(i,j,m)+ud(i,j,m)*ud(i,j,m))
                                    !Angle of (1,0) geospace vector 
                                    !in Mercator projection
        v2=ATAN2(dy(i,j),dx(i,j))
                                    !Angle in Mercator space
        v2=v2+v1
        ud(i,j,m)=dl*COS(v2) !u in Mercator
        vd(i,j,m)=dl*SIN(v2) !v in Mercator
       ENDDO
      ENDDO
      write(*,*)'ud,vd rotated from geospace to Mercator'
      WRITE(1,100)((ud(i,j,1),i=1,nx),j=1,ny) 
      WRITE(1,100)((vd(i,j,1),i=1,nx),j=1,ny) 
      ENDIF !geo_to_merc
c
cc
cc      stop
c --- ------------------------------------------------------------
      STOP
      END    


          

