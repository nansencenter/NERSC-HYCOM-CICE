      SUBROUTINE rotate(ud,vd,mlat,mlon,nx,ny)
c ---------------------------------------------------------------------
c --- Author: Knud Simonsen, NERSC, 01.04.1996 
c --- Modified:
c ---------------------------------------------------------------------
c --- Rotates a vector field from a grid in geographical coordinates
c --- into the grid defined by the lat,lon in the input variables.
c --- C-grid is assumed.                     
c ---  
c ---  Coordinate  system                   | i
c ---                                       |
c ---                                       |
c ---                                       |
c ---                           j-----------|
c ---------------------------------------------------------------------
c --- Input: mlat, mlon: position in scalar point
c ---        nx,ny:     dimension of the model grid.
c ---        ud,vd:  Unrotated vector components,  where ud is the EW
c ---                component and vd is the NS component 
c --- Output: ud,vd: Rotated vector components, where  ud is along the
c ---                i-axis and vd is along the j-axis.
c ----------------------------------------------------------------------
      INTEGER nx,ny,i,j,mo,im,jm,ip,jp
      REAL ud(nx,ny,12),vd(nx,ny,12)
     +     ,mlat(nx,ny),mlon(nx,ny)
     +     ,pi,pi2,radian,radinv,fli,flj
     +     ,u_up,v_up,u_vp,v_vp,theta_up,theta_vp
     +     ,dlon,dd(0:122,0:72),dlat

      data radian/57.29578/,pi/3.14159265/
      pi2 = pi/2.
      radinv=1./radian

c --- ---------------------------------------- test wind -------
c      DO i=1,nx
c       DO j=1,ny
c        ud(i,j,1)=3. 
c        vd(i,j,1)=4. 
c       ENDDO
c      ENDDO
ccc
c      write(*,*)'ud,vd test values'
c      CALL dmp2dtec(ud,iflg,dd,nx,ny,'r')
c      CALL dmp2dtec(vd,iflg,dd,nx,ny,'r')
cc      WRITE(1,100)((ud(i,j,1),i=1,nx),j=1,ny) 
c      WRITE(1,100)((vd(i,j,1),i=1,nx),j=1,ny) 
cc.test      GOTO 900
      WRITE(*,*)'------------------------------------------------'
      WRITE(*,*)'Rotate wind vectores and interpolate into C-grid'
      WRITE(*,*)'(Subroutine rotate.f)'
      WRITE(*,*)'------------------------------------------------'

c --- ----------------------------------------------------------
c --- Assumes that all parameters are provided in scalar point 
c --- and interpolates into the U- and V (C-grid) points, and
c --- perform the rotation rquired in curvlinear grid.
c --- ----------------------------------------------------------                             
      DO j=ny,1,-1
       DO i=nx,1,-1
        im=MAX(i-1,1)
        fli=FLOAT(i-im)
        jm=MAX(j-1,1)
        flj=FLOAT(j-jm)
                                       !Rotation angle in u-point 
        dlon=(mlon(i,j)-mlon(im,j))
        dlat=(mlat(i,j)-mlat(im,j))
        IF(dlon.LT.180.)dlon=360.0+dlon
        IF(dlon.GT.180.)dlon=dlon-360.0
        theta_up = atan2(dlat, 
     &            dlon*cos(radinv*.5*(mlat(i,j)+mlat(im,j))) )

                                       !Rotation angle in v-point 
        dlon=(mlon(i,j)-mlon(i,jm))
        dlat=mlat(i,j)-mlat(i,jm)

        IF(dlon.LT.180.)dlon=360.0+dlon
        IF(dlon.GT.180.)dlon=dlon-360.0
        theta_vp = atan2(dlat,
     &            dlon*cos(radinv*.5*(mlat(i,j)+mlat(i,jm))) )

        DO mo=1,12 
         u_up=.5*(ud(i,j,mo)+ud(im,j,mo)) !Unrotated vel. in u-point
         v_up=.5*(vd(i,j,mo)+vd(im,j,mo))
 
         u_vp=.5*(ud(i,j,mo)+ud(i,jm,mo)) !Unrotated vel. in v-point
         v_vp=-.5*(vd(i,j,mo)+vd(i,jm,mo))

                        
                                          !Perform rotation:
                                          !In U-point only
c         ud(i,j,mo)= (u_up*COS(theta_up)+ v_up*SIN(theta_up))*fli
c         vd(i,j,mo)=(-u_up*SIN(theta_up)+ v_up*COS(theta_up))*flj

                                          !In V-point only

c         ud(i,j,mo)= (u_vp*SIN(theta_vp)+ v_vp*COS(theta_vp))*fli
c         vd(i,j,mo)=-(-u_vp*COS(theta_vp)+ v_vp*SIN(theta_vp))*flj

                                          !In C-grid
         ud(i,j,mo)= (u_up*COS(theta_up)+ v_up*SIN(theta_up))*fli
         vd(i,j,mo)=-(-u_vp*COS(theta_vp)+ v_vp*SIN(theta_vp))*flj

        ENDDO  !i-loop
       ENDDO  !j-loop
      END DO  !mo-loop
c      WRITE(1,100)((ud(i,j,1),i=1,nx),j=1,ny)
c      WRITE(1,100)((vd(i,j,1),i=1,nx),j=1,ny)
c      CALL dmp2dtec(ud,iflg,dd,nx,ny,'r')
c      CALL dmp2dtec(vd,iflg,dd,nx,ny,'r')
   99 FORMAT(20I4)
  100 FORMAT(10(1x,e12.6))
      RETURN
      END

