      SUBROUTINE curlin_geo(ud,vd,mlon,mlat,ulon,ulat
     &                     ,vlon,vlat,nx,ny,nz)
c ---------------------------------------------------------------------
c --- Author: Knud Simonsen, NERSC, 26.06.1996
c --- Modified:
c --- ------------------------------------------------------------------------   
c --- Assumes that ud,vd contains velocities in xi, eta space on C-grid
c --- and rotates to geospace. 
c --- ------------------------------------------------------------------------   
      INTEGER nx,ny,nz,i,j,mo,im,jm,ip,jp
      REAL ud(nx,ny,nz),vd(nx,ny,nz)      
     +     ,ulat(nx,0:ny),ulon(nx,0:ny)
     +     ,vlat(0:nx,ny),vlon(0:nx,ny)
     +     ,mlat(nx,ny),mlon(nx,ny)
     +     ,u_up,v_up,u_vp,v_vp,theta_up,theta_vp
     +     ,dlon,v1,v2,dl,u_pp,v_pp
     +     ,dr,dra(4)

      LOGICAL polcorr
 
      data radian/57.29578/,pi/3.14159265/
      pi2 = pi/2.
      radinv=1./radian

      polcorr = .False.          !Correction around the Pole.   





      DO j=1,ny
       DO i=1,nx
        ip=MIN(i+1,nx)
        fli=FLOAT(i-ip)
        jp=MIN(j+1,ny)
        flj=FLOAT(j-jp)


                                       !Rotation of u in p-point
        dlon=ulon(ip,j)-ulon(i,j)
        IF(dlon.LT.-180.)dlon=360.0+dlon
        IF(dlon.GT.180.)dlon=dlon-360.0
        theta_up = atan2(ulat(ip,j)-ulat(i,j),
     &            dlon*cos(radinv*.5*(ulat(ip,j)+ulat(i,j))) )

                                       !Rotation of v in p-point
        dlon=vlon(i,jp)-vlon(i,j)
        IF(dlon.LT.-180.)dlon=360.0+dlon
        IF(dlon.GT.180.)dlon=dlon-360.0
        theta_vp = atan2(vlat(i,jp)-vlat(i,j),
     &            dlon*cos(radinv*.5*(vlat(i,jp)+vlat(i,j))) )
     &            -pi2 
c.diag        IF(ABS(ulon(ip,j)).GT.170)
c.diag     &  write(*,*)i,j,lon,dlon,theta_vp,theta_up
         
        DO mo=1,nz
c         u_pp=.5*(ud(i,j,mo)+ud(ip,j,mo))  !u in p-point
c         v_pp=.5*(vd(i,j,mo)+vd(i,jp,mo))  !v in p-point
                                            !ud,vd given in p-point
          u_pp= ud(i,j,mo)
          v_pp= vd(i,j,mo)

                                           !Rotate to geo. space
         ud(i,j,mo)=u_pp*COS(theta_up)-v_pp*SIN(theta_vp)
         vd(i,j,mo)=u_pp*SIN(theta_up)+v_pp*COS(theta_vp)
c         ud(i,j,mo)=theta_up
c         vd(i,j,mo)=theta_vp 
        ENDDO
       ENDDO
      ENDDO

      IF(polcorr)THEN
                                          !Perform interpolation over the
                                          !polar area in order to hide the
                                          !error arising from the calc. of
                                          !the rot. angle
      mo=1
      im=103
      ip=109
      jm=42
      jp=55 
      write(*,*)'CoRrection at the Pole' 
      DO i=im+1,ip-1
       DO j=jm+1,jp-1
        ud(i,j,mo)=0.
        vd(i,j,mo)=0.
        dum=0.
        DO ii=im,ip
         DO jj=jm,jp,jp-jm
          dr=111200./spherdist(mlon(ii,jj),mlat(ii,jj)
     &                        ,mlon(i,j),mlat(i,j))
          dr=ABS(dr)
c          write(*,*)ii,jj,i,j,dr
          ud(i,j,mo)=ud(i,j,mo)+dr*ud(ii,jj,mo)
          vd(i,j,mo)=vd(i,j,mo)+dr*vd(ii,jj,mo)
          dum=dum+dr
         ENDDO
        ENDDO 
        DO ii=im,ip,ip-im
         DO jj=jm+1,jp-1,1
          dr=111200./spherdist(mlon(ii,jj)*radian,mlat(ii,jj)*radian
     &                    ,mlon(i,j)*radian,mlat(i,j)*radian)
          dr=ABS(dr)
c          write(*,*)ii,jj,i,j,dr
          ud(i,j,mo)=ud(i,j,mo)+dr*ud(ii,jj,mo)
          vd(i,j,mo)=vd(i,j,mo)+dr*vd(ii,jj,mo)
          dum=dum+dr
         ENDDO
        ENDDO 
        ud(i,j,mo)=ud(i,j,mo)/dum
        vd(i,j,mo)=vd(i,j,mo)/dum
       ENDDO
      ENDDO
      ENDIF

      RETURN
      END



 
