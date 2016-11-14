SUBROUTINE rotate2(ud,vd,mlat,mlon,dir,keeppoint)
! --- ----------------------------------------------------------------------
! --- Rotates a vector field from a grid in geographical coordinates
! --- into or from the grid defined by the lat,lon in the input variables.
! --- C-grid is assumed. If keeppoint is true, in/out velocities are 
! --- given in the same point. 
! ---  
! --- Input: mlat, mlon: position in scalar point
! ---        nx,ny     : dimension of the model grid.
! ---        ud,vd     : Unrotated vector components,  where ud is the EW
! ---                    component and vd is the NS component 
! ---        dir       : l2m (latlon to general)
! ---                    m2l (general to latlon)
! ---        keeppoint : keep vector in same point!)
! --- 
! --- Output: ud,vd: Rotated vector components, where  ud is along the
! ---               i-axis and vd is along the j-axis.
! --- ----------------------------------------------------------------------
    use mod_xc
      implicit none
      real, intent(inout), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: ud,vd
      real, intent(in   ), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: mlon, mlat
      character(len=3), intent(in) :: dir
      logical         , intent(in) :: keeppoint
      integer i,j
      real pi,pi2,radian,radinv
      real u_up,v_up,u_vp,v_vp,theta_up,theta_vp,up,vp
      real dlon,dlat
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: urot, vrot
      data radian/57.29578/,pi/3.14159265/
      pi2 = pi/2.
      radinv=1./radian
!
! ---  ----------------------------------------------------------
! ---  Assumes that all parameters are provided in scalar point 
! ---  and interpolates into the U- and V (C-grid) points, and
! ---  perform the rotation rquired in curvlinear grid.
! ---  -------------------------------------------------------   
!
! --- Goes from u/v in scalar points to u/v in velocity points
      urot=0.0
      vrot=0.0
      if (dir == 'l2m') then
         do j=1-nbdy+1,jdm+nbdy-1
         do i=1-nbdy+1,idm+nbdy-1
! TODO      Correct rotation for keeppoint
! ---       Rotation angle in u-point 
            dlon=(mlon(i,j)-mlon(i-1,j))
            dlat=(mlat(i,j)-mlat(i-1,j))
            if(dlon.LT.180.)dlon=360.0+dlon
            if(dlon.GT.180.)dlon=dlon-360.0
            theta_up = atan2(dlat,dlon*cos(radinv*.5*(mlat(i,j)+mlat(i-1,j))) )
! 
! ---       Rotation angle in v-point 
            dlon=mlon(i,j)-mlon(i,j-1)
            dlat=mlat(i,j)-mlat(i,j-1)
!
            if(dlon.LT.180.)dlon=360.0+dlon
            if(dlon.GT.180.)dlon=dlon-360.0
            theta_vp = atan2(dlat,dlon*cos(radinv*.5*(mlat(i,j)+mlat(i,j-1))) )
!
            if (keeppoint) then
! ---          Unrotated vel. in original point
               u_up= ud(i,j)
               u_vp= ud(i,j)
               v_up= vd(i,j)
               v_vp= vd(i,j)
            else
! ---          Unrotated vel. in u-point
               u_up=.5*(ud(i,j)+ud(i-1,j))
               v_up=.5*(vd(i,j)+vd(i-1,j))
! ---          Unrotated vel. in v-point
               u_vp=.5*(ud(i,j)+ud(i,j-1))
               v_vp=.5*(vd(i,j)+vd(i,j-1))
            end if
!
! ---       Final rotated velocities
            urot(i,j)= u_up*COS(theta_up)+ v_up*SIN(theta_up)
            vrot(i,j)= u_vp*COS(theta_vp)+ v_vp*SIN(theta_vp)
         enddo

         enddo
!
! --- Goes from u/v in velocity points to u/v in scalar points
      elseif (dir == 'm2l') then
         do j=1-nbdy+1,jdm+nbdy-1
         do i=1-nbdy+1,idm+nbdy-1
!
! ---       Rotation angle in p-point 
            dlon=mlon(i+1,j)-mlon(i-1,j)
            dlat=mlat(i+1,j)-mlat(i-1,j)
            if(dlon.LT.180.)dlon=360.0+dlon
            if(dlon.GT.180.)dlon=dlon-360.0
            theta_up = atan2(dlat,dlon*cos(radinv*.5*(mlat(i-1,j)+mlat(i+1,j))) )
!    
! ---       Rotation angle in p-point 
            dlon=mlon(i,j+1)-mlon(i,j-1)
            dlat=mlat(i,j+1)-mlat(i,j-1)
            if(dlon.LT.180.)dlon=360.0+dlon
            if(dlon.GT.180.)dlon=dlon-360.0
            theta_vp = atan2(dlat,dlon*cos(radinv*.5*(mlat(i,j-1)+mlat(i,j+1))) )
! 
            if (keeppoint) then
! ---          Unrotated vel. in original point
               up=ud(i,j)
               vp=vd(i,j)
            else
! ---          Unrotated vel. in p-point
               up=.5*(ud(i,j)+ud(i+1,j))
               vp=.5*(vd(i,j)+vd(i,j+1))
            end if
!
! ---       Final rotated velocities
            urot(i,j)= up*cos(theta_up)+ vp*cos(theta_vp)
            vrot(i,j)= up*sin(theta_up)+ vp*sin(theta_vp)
         enddo
         enddo
      else
         if (mnproc==1) then
            write(lp,'(a)') 'Unknown rotation dir '//dir
         end if
         call xcstop('(mod_hycom_nersc:rotate)')
         stop '(mod_hycom_nersc:rotate)'
      endif
      ud=urot
      vd=vrot
END subroutine rotate2
