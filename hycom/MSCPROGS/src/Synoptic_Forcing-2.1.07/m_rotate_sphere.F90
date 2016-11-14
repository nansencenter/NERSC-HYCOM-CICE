module m_rotate_sphere

private :: cross_product, geo2cart,norm2 
contains
SUBROUTINE rotate_sphere(ud,vd,mlat,mlon,nx,ny,dir)
! Rotates a vector field from a grid in geographical coordinates
! into or from the grid defined by the lat,lon in the input variables.
! C-grid is assumed.                     
!  
! Input: mlat, mlon: position in scalar point
!        nx,ny:     dimension of the model grid.
!        ud,vd:  Unrotated vector components,  where ud is the EW
!                component and vd is the NS component 
!        dir:    l2m (latlon to general)
!                m2l (general to latlon)
!
! Output: ud,vd: Rotated vector components, where  ud is along the
!                i-axis and vd is along the j-axis.
!
! KAL - This version shouldbe more accurate near the geo poles
! ----------------------------------------------------------------------
   implicit none
   integer, intent(in) :: nx
   integer, intent(in) :: ny
   real, intent(inout) :: ud(nx,ny),vd(nx,ny)
   real, intent(in)    :: mlon(nx,ny)
   real, intent(in)    :: mlat(nx,ny)
   character(len=3), intent(in) :: dir

   integer i,j,im,jm,ip,jp
   real   u_up,v_up,u_vp,v_vp,theta_up,theta_vp,up,vp
   real   dlon,dlat
   real :: urot(nx,ny),vrot(nx,ny)
   real*8 pi,radian,radinv
   real*8, dimension(3) :: a,b,c,n1,n2,n3,un3
   real*8 :: absang, signang
   real   :: blon,blat

   data radian/57.29578/,pi/3.14159265/
   radinv=1./radian


! ----------------------------------------------------------
! Assumes that all parameters are provided in scalar point 
! and interpolates into the U- and V (C-grid) points, and
! perform the rotation rquired in curvlinear grid.
! -------------------------------------------------------   
   

   urot=0.0 
   vrot=0.0
   if (dir == 'l2m') then
!     !$OMP PARALLEL DO PRIVATE (i,j,dlon,dlat,theta_up,theta_vp,u_up,v_up, &
!     !$OMP                      u_vp,v_vp,up,vp) &
!     !$OMP SCHEDULE(STATIC,jblk)
      do j=2,ny
         do i=2,nx
! Rotation angle in u-point 
            !dlon=(mlon(i,j)-mlon(i-1,j))
            !dlat=(mlat(i,j)-mlat(i-1,j))
            !if(dlon.LT.180.)dlon=360.0+dlon
            !if(dlon.GT.180.)dlon=dlon-360.0
            !theta_up = atan2(dlat,dlon*cos(radinv*.5*(mlat(i,j)+mlat(i-1,j))) )


            ! This is used to  defines the great circle which is tangent to the ! longitude line at this point
            blon=mlon(i,j)+90. ; blat=0.  ! NB - check sign of this

            ! 3D position of points
            a=geo2cart(mlon(i,j),mlat(i,j))
            b=geo2cart(blon,blat)
            c=geo2cart(mlon(i+1,j),mlat(i+1,j))

            ! normal to the plane spanned by a and b
            n1=cross_product(a,b)
            n2=cross_product(a,c)

            ! normal to the normals
            n3=cross_product(n1,n2)

            ! unit normal
            un3=n3/norm2(n3,3)

            ! Sign of normal vector rel NP vector
            signang=sign(1._8,dot_product(geo2cart(0.,90.),un3))

            ! Change orientation on S hemisphere
            signang=signang*sign(1.,mlat(i,j));


            theta_up=atan2(signang*dot_product(n3,un3),dot_product(n1,n2))

    
! Rotation angle in v-point 

            !dlon=(mlon(i,j)-mlon(i,j-1))
            !dlat=mlat(i,j)-mlat(i,j-1)
            !if(dlon.LT.180.)dlon=360.0+dlon
            !if(dlon.GT.180.)dlon=dlon-360.0
            !theta_vp = atan2(dlat,dlon*cos(radinv*.5*(mlat(i,j)+mlat(i,j-1))) )

            ! This is used to  defines the great circle which is tangent to the ! longitude line at this point
            blon=mlon(i,j)+90. ; blat=0.  ! NB - check sign of this

            ! 3D position of points
            a=geo2cart(mlon(i,j),mlat(i,j))
            b=geo2cart(blon,blat)
            c=geo2cart(mlon(i,j+1),mlat(i,j+1))

            ! normal to the plane spanned by a and b
            n1=cross_product(a,b)
            n2=cross_product(a,c)

            ! normal to the normals
            n3=cross_product(n1,n2)

            ! unit normal
            un3=n3/norm2(n3,3)

            ! Sign of normal vector rel sp vector
            signang=sign(1._8,dot_product(geo2cart(0.,90.),un3))

            ! Change orientation on S hemisphere
            signang=signang*sign(1.,mlat(i,j));

            theta_vp=atan2(signang*dot_product(n3,un3),dot_product(n1,n2))

            !print *,i,j,theta_up,theta_vp,signang

! Final calculation of rotated vectors


            u_up=.5*(ud(i,j)+ud(i-1,j)) !Unrotated vel. in u-point
            v_up=.5*(vd(i,j)+vd(i-1,j))
    
            u_vp=.5*(ud(i,j)+ud(i,j-1)) !Unrotated vel. in v-point
            v_vp=.5*(vd(i,j)+vd(i,j-1))

!           u_up=(ud(i,j)) !Unrotated vel. in u-point
!           v_up=(vd(i,j))
!   
!           u_vp=(ud(i,j)) !Unrotated vel. in v-point
!           v_vp=(vd(i,j))
!
!            print *,'4',theta_up,theta_vp,u_up,u_vp,v_up,v_vp

            urot(i,j)= u_up*COS(theta_up)+ v_up*SIN(theta_up)
            vrot(i,j)= u_vp*COS(theta_vp)+ v_vp*SIN(theta_vp)
         enddo 
      enddo
      !$OMP END PARALLEL DO


   elseif (dir == 'm2l') then
      stop '(m2l step not finished for this version of rotate)'
      !$OMP PARALLEL DO PRIVATE (i,j,dlon,dlat,theta_up,theta_vp,u_up,v_up, &
      !$OMP                      u_vp,v_vp,up,vp) &
      !$OMP SCHEDULE(STATIC,jblk)
      do j=2,ny-1
         do i=2,nx-1
! Rotation angle in p-point 
            dlon=mlon(i+1,j)-mlon(i-1,j)
            dlat=mlat(i+1,j)-mlat(i-1,j)
            if(dlon.LT.180.)dlon=360.0+dlon
            if(dlon.GT.180.)dlon=dlon-360.0
            theta_up = atan2(dlat,dlon*cos(radinv*.5*(mlat(i-1,j)+mlat(i+1,j))) )
    
! Rotation angle in p-point 
            dlon=mlon(i,j+1)-mlon(i,j-1)
            dlat=mlat(i,j+1)-mlat(i,j-1)
            if(dlon.LT.180.)dlon=360.0+dlon
            if(dlon.GT.180.)dlon=dlon-360.0
            theta_vp = atan2(dlat,dlon*cos(radinv*.5*(mlat(i,j-1)+mlat(i,j+1))) )

!Unrotated vel. in p-point
            up=.5*(ud(i,j)+ud(i+1,j))
            vp=.5*(vd(i,j)+vd(i,j+1))
    
            urot(i,j)= up*cos(theta_up)+ vp*cos(theta_vp)
            vrot(i,j)= up*sin(theta_up)+ vp*sin(theta_vp)
         enddo 
      enddo
      !$OMP END PARALLEL DO
   else
      stop 'rotate_sphere'
   endif
   ud=urot
   vd=vrot

END subroutine rotate_sphere

   function geo2cart(lon,lat)
      implicit none
      real*8, parameter :: rad=1.7453292519943295E-02,deg=57.29577951308232
      real  , intent(in) :: lon,lat
      real*8, dimension(3) :: geo2cart
      ! Assume radius of earth == 1

      real*8 :: lambda,theta

      lambda=lat*rad
      theta=lon*rad
      geo2cart(1)=cos(lambda)*cos(theta)
      geo2cart(2)=cos(lambda)*sin(theta)
      geo2cart(3)=sin(lambda)
   end function geo2cart


   ! Routine to calculate cross product of two 3D vectors.
   function cross_product(v1,v2)
      implicit none
      real*8, intent(in), dimension(3) :: v1,v2
      real*8, dimension(3) :: cross_product

      cross_product(1) = v1(2)*v2(3) - v1(3)*v2(2)
      cross_product(2) = v1(3)*v2(1) - v1(1)*v2(3)
      cross_product(3) = v1(1)*v2(2) - v1(2)*v2(1)
   end function cross_product

   ! Routine to calculate vector norm (more precise the 2-norm)
   function norm2(vector,p)
      implicit none
      integer, intent(in) :: p
      real*8,    intent(in) :: vector(p)
      integer :: i
      real*8    :: norm2

      norm2=0.
      do i=1,p
         norm2 = norm2 + vector(i)**2
      end do
      norm2=sqrt(norm2)
   end function norm2



end module m_rotate_sphere
