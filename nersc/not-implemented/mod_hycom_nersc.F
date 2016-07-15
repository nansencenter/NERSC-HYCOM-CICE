      module mod_hycom_nersc
c --- Module contains various routines with no dependencies 
c --- to other modules  than mod_xc and common_blocks
c --- These routines are called by various nersc processing routines
      implicit none
c --- Some convenient variables to keep here
      integer         ,save :: imem
      character(len=3),save :: rungen   ! version code for run
      logical,         save :: ensflag  ! Ensemble run flag
      private :: ncerr, ncdraft_1
      contains



c --- Print status to screen and to output file $rungenstatus
      subroutine prtstatus(dtime)
         use mod_xc
         implicit none
         real*8, external :: wtime
#if ! defined (MPI) && ! defined (SHMEM)
         integer, external :: omp_get_max_threads
#endif
         real*8,  intent(in) :: dtime
         real*8, save :: time_last=0.0,  time_now=0.0,
     &      tot_time=0.0, time_per_tstep=0.0
         logical, save :: first=.true.
         integer, save :: steps=0, icpu = 1
         integer :: iyear, iday, ihour, iss
         include 'common_blocks.h'
c
         if (first) then
            first=.false.
#if ! defined (MPI) && ! defined (SHMEM)
            icpu = max(1,omp_get_max_threads())
#else
            icpu=1
#endif
            time_now=wtime()
            time_last=time_now
         else
            time_last=time_now
            time_now=wtime()
         endif
         steps=steps+1
         tot_time = tot_time + (time_now-time_last)
         time_per_tstep=tot_time/steps
c
         call forday(dtime,yrflag,iyear,iday,ihour)
         iss=nint((dtime-floor(dtime,kind=8))*86400.d0) - ihour*3600
c
         if (mnproc==1) then
           open(10,file=rungen//'status',status='unknown')
           write(10,1000) imem,nstep,iyear,iday-1,ihour,iss,
     &        time_now-time_last,time_per_tstep,tot_time/3600.0d0
           close(10)
           write(lp,1000) imem,nstep,iyear,iday-1,ihour,iss,
     &        time_now-time_last,time_per_tstep,tot_time/3600.0d0
           call flush(lp)
         end if
1000  format('iens=',i4,' nstep=',i10,' --- y',i4.4,'_d',i3.3,'_h',i2.2,
     &       '_s',i4.4,' --- CPU=',f6.2,', Avg CPU=',f6.2,' Tot CPU=',
     &       f6.2)
      end subroutine prtstatus
c
c
c
c --- ----------------------------------------------------------------------
c --- Rotates a vector field from a grid in geographical coordinates
c --- into or from the grid defined by the lat,lon in the input variables.
c --- C-grid is assumed.                     
c ---  
c --- Input: mlat, mlon: position in scalar point
c ---        nx,ny     : dimension of the model grid.
c ---        ud,vd     : Unrotated vector components,  where ud is the EW
c ---                    component and vd is the NS component 
c ---        dir       : l2m (latlon to general)
c ---                    m2l (general to latlon)
c --- 
c --- Output: ud,vd: Rotated vector components, where  ud is along the
c ---               i-axis and vd is along the j-axis.
c --- ----------------------------------------------------------------------
      SUBROUTINE rotate(ud,vd,mlat,mlon,nx,ny,dir)
      use mod_xc
      implicit none
      integer, intent(in) :: nx
      integer, intent(in) :: ny
      real, intent(inout) :: ud(nx,ny),vd(nx,ny)
      real, intent(in)    :: mlon(nx,ny)
      real, intent(in)    :: mlat(nx,ny)
      character(len=3), intent(in) :: dir
c
      integer i,j
      real pi,pi2,radian,radinv
      real u_up,v_up,u_vp,v_vp,theta_up,theta_vp,up,vp
      real dlon,dlat
      real :: urot(nx,ny),vrot(nx,ny)
      data radian/57.29578/,pi/3.14159265/
      pi2 = pi/2.
      radinv=1./radian
c
c ---  ----------------------------------------------------------
c ---  Assumes that all parameters are provided in scalar point 
c ---  and interpolates into the U- and V (C-grid) points, and
c ---  perform the rotation rquired in curvlinear grid.
c ---  -------------------------------------------------------   
c
      urot=0.0 
      vrot=0.0
      if (dir == 'l2m') then
!$OMP PARALLEL DO PRIVATE (i,j,dlon,dlat,theta_up,theta_vp,u_up,v_up, 
!$OMP&                     u_vp,v_vp,up,vp) 
!$OMP&SCHEDULE(STATIC,jblk)
         do j=2,ny
         do i=2,nx
c ---       Rotation angle in u-point 
            dlon=(mlon(i,j)-mlon(i-1,j))
            dlat=(mlat(i,j)-mlat(i-1,j))
            if(dlon.LT.180.)dlon=360.0+dlon
            if(dlon.GT.180.)dlon=dlon-360.0
            theta_up = atan2(dlat,
     &         dlon*cos(radinv*.5*(mlat(i,j)+mlat(i-1,j))) )
c 
c ---       Rotation angle in v-point 
            dlon=(mlon(i,j)-mlon(i,j-1))
            dlat=mlat(i,j)-mlat(i,j-1)
c
            if(dlon.LT.180.)dlon=360.0+dlon
            if(dlon.GT.180.)dlon=dlon-360.0
            theta_vp = atan2(dlat,
     &         dlon*cos(radinv*.5*(mlat(i,j)+mlat(i,j-1))) )
c
c ---       Unrotated vel. in u-point
            u_up=.5*(ud(i,j)+ud(i-1,j))
            v_up=.5*(vd(i,j)+vd(i-1,j))
c ---       Unrotated vel. in v-point
            u_vp=.5*(ud(i,j)+ud(i,j-1)) 
            v_vp=.5*(vd(i,j)+vd(i,j-1))
c
c ---       Final rotated velocities
            urot(i,j)= u_up*COS(theta_up)+ v_up*SIN(theta_up)
            vrot(i,j)= u_vp*COS(theta_vp)+ v_vp*SIN(theta_vp)
         enddo 
         enddo
C$OMP END PARALLEL DO
c
      elseif (dir == 'm2l') then
c$OMP PARALLEL DO PRIVATE (i,j,dlon,dlat,theta_up,theta_vp,u_up,v_up, 
c$OMP&                     u_vp,v_vp,up,vp) 
c$OMP&SCHEDULE(STATIC,jblk)
         do j=2,ny-1
         do i=2,nx-1
c
c ---       Rotation angle in p-point 
            dlon=mlon(i+1,j)-mlon(i-1,j)
            dlat=mlat(i+1,j)-mlat(i-1,j)
            if(dlon.LT.180.)dlon=360.0+dlon
            if(dlon.GT.180.)dlon=dlon-360.0
            theta_up = atan2(dlat,
     &         dlon*cos(radinv*.5*(mlat(i-1,j)+mlat(i+1,j))) )
c    
c ---       Rotation angle in p-point 
            dlon=mlon(i,j+1)-mlon(i,j-1)
            dlat=mlat(i,j+1)-mlat(i,j-1)
            if(dlon.LT.180.)dlon=360.0+dlon
            if(dlon.GT.180.)dlon=dlon-360.0
            theta_vp = atan2(dlat,
     &         dlon*cos(radinv*.5*(mlat(i,j-1)+mlat(i,j+1))) )
c 
c ---       Unrotated vel. in p-point
            up=.5*(ud(i,j)+ud(i+1,j))
            vp=.5*(vd(i,j)+vd(i,j+1))
c
c ---       Final rotated velocities
            urot(i,j)= up*cos(theta_up)+ vp*cos(theta_vp)
            vrot(i,j)= up*sin(theta_up)+ vp*sin(theta_vp)
         enddo 
         enddo
c$OMP END PARALLEL DO
      else
         if (mnproc==1) then
            write(lp,'(a)') 'Unknown rotation dir '//dir
         end if
         call xcstop('(mod_hycom_nersc:rotate)')
         stop '(mod_hycom_nersc:rotate)'
      endif
      ud=urot
      vd=vrot
      END subroutine rotate

c --- ----------------------------------------------------------------------
c --- Rotates a vector field from a grid in geographical coordinates
c --- into or from the grid defined by the lat,lon in the input variables.
c --- C-grid is assumed. If keeppoint is true, in/out velocities are 
c --- given in the same point. 
c ---  
c --- Input: mlat, mlon: position in scalar point
c ---        nx,ny     : dimension of the model grid.
c ---        ud,vd     : Unrotated vector components,  where ud is the EW
c ---                    component and vd is the NS component 
c ---        dir       : l2m (latlon to general)
c ---                    m2l (general to latlon)
c ---        keeppoint : keep vector in same point!)
c --- 
c --- Output: ud,vd: Rotated vector components, where  ud is along the
c ---               i-axis and vd is along the j-axis.
c --- ----------------------------------------------------------------------
      SUBROUTINE rotate2(ud,vd,mlat,mlon,dir,keeppoint)
      use mod_xc
      implicit none
      real, intent(inout), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &   ud,vd
      real, intent(in   ), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &   mlon, mlat
      character(len=3), intent(in) :: dir
      logical         , intent(in) :: keeppoint
c
      integer i,j
      real pi,pi2,radian,radinv
      real u_up,v_up,u_vp,v_vp,theta_up,theta_vp,up,vp
      real dlon,dlat
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: urot, vrot
      data radian/57.29578/,pi/3.14159265/
      pi2 = pi/2.
      radinv=1./radian
c
c ---  ----------------------------------------------------------
c ---  Assumes that all parameters are provided in scalar point 
c ---  and interpolates into the U- and V (C-grid) points, and
c ---  perform the rotation rquired in curvlinear grid.
c ---  -------------------------------------------------------   
c
c --- Goes from u/v in scalar points to u/v in velocity points
      urot=0.0 
      vrot=0.0
      if (dir == 'l2m') then
!$OMP PARALLEL DO PRIVATE (i,j,dlon,dlat,theta_up,theta_vp,u_up,v_up, 
!$OMP&                     u_vp,v_vp,up,vp) 
!$OMP&SCHEDULE(STATIC,jblk)
         do j=1-nbdy+1,jj+nbdy-1
         do i=1-nbdy+1,ii+nbdy-1
c TODO      Correct rotation for keeppoint
c ---       Rotation angle in u-point 
            dlon=(mlon(i,j)-mlon(i-1,j))
            dlat=(mlat(i,j)-mlat(i-1,j))
            if(dlon.LT.180.)dlon=360.0+dlon
            if(dlon.GT.180.)dlon=dlon-360.0
            theta_up = atan2(dlat,
     &         dlon*cos(radinv*.5*(mlat(i,j)+mlat(i-1,j))) )
c 
c ---       Rotation angle in v-point 
            dlon=mlon(i,j)-mlon(i,j-1)
            dlat=mlat(i,j)-mlat(i,j-1)
c
            if(dlon.LT.180.)dlon=360.0+dlon
            if(dlon.GT.180.)dlon=dlon-360.0
            theta_vp = atan2(dlat,
     &         dlon*cos(radinv*.5*(mlat(i,j)+mlat(i,j-1))) )
c
            if (keeppoint) then
c ---          Unrotated vel. in original point
               u_up= ud(i,j)
               u_vp= ud(i,j)
               v_up= vd(i,j)
               v_vp= vd(i,j)
            else
c ---          Unrotated vel. in u-point
               u_up=.5*(ud(i,j)+ud(i-1,j))
               v_up=.5*(vd(i,j)+vd(i-1,j))
c ---          Unrotated vel. in v-point
               u_vp=.5*(ud(i,j)+ud(i,j-1)) 
               v_vp=.5*(vd(i,j)+vd(i,j-1))
            end if
c
c ---       Final rotated velocities
            urot(i,j)= u_up*COS(theta_up)+ v_up*SIN(theta_up)
            vrot(i,j)= u_vp*COS(theta_vp)+ v_vp*SIN(theta_vp)
         enddo 
         enddo
C$OMP END PARALLEL DO
c
c --- Goes from u/v in velocity points to u/v in scalar points
      elseif (dir == 'm2l') then
c$OMP PARALLEL DO PRIVATE (i,j,dlon,dlat,theta_up,theta_vp,u_up,v_up, 
c$OMP&                     u_vp,v_vp,up,vp) 
c$OMP&SCHEDULE(STATIC,jblk)
         do j=1-nbdy+1,jj+nbdy-1
         do i=1-nbdy+1,ii+nbdy-1
c
c ---       Rotation angle in p-point 
            dlon=mlon(i+1,j)-mlon(i-1,j)
            dlat=mlat(i+1,j)-mlat(i-1,j)
            if(dlon.LT.180.)dlon=360.0+dlon
            if(dlon.GT.180.)dlon=dlon-360.0
            theta_up = atan2(dlat,
     &         dlon*cos(radinv*.5*(mlat(i-1,j)+mlat(i+1,j))) )
c    
c ---       Rotation angle in p-point 
            dlon=mlon(i,j+1)-mlon(i,j-1)
            dlat=mlat(i,j+1)-mlat(i,j-1)
            if(dlon.LT.180.)dlon=360.0+dlon
            if(dlon.GT.180.)dlon=dlon-360.0
            theta_vp = atan2(dlat,
     &         dlon*cos(radinv*.5*(mlat(i,j-1)+mlat(i,j+1))) )
c 
            if (keeppoint) then
c ---          Unrotated vel. in original point
               up=ud(i,j)
               vp=vd(i,j)
            else
c ---          Unrotated vel. in p-point
               up=.5*(ud(i,j)+ud(i+1,j))
               vp=.5*(vd(i,j)+vd(i,j+1))
            end if
c
c ---       Final rotated velocities
            urot(i,j)= up*cos(theta_up)+ vp*cos(theta_vp)
            vrot(i,j)= up*sin(theta_up)+ vp*sin(theta_vp)
         enddo 
         enddo
c$OMP END PARALLEL DO
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











c
c
c
C --- -----------------------------------------
C --- Computes the distance between geo. pos.
C --- lon1,lat1 and lon2,lat2. 
C --- INPUT is in degrees.
C --- TODO: Replace with Haversine formula (more efficient
C --- and better precision compared to acos operation)
C --- -----------------------------------------
      real function spherdist(lon1,lat1,lon2,lat2)
      implicit none
      REAL, intent(in) :: lon1,lat1,lon2,lat2 ! Pos. in degrees
c
      real*8, parameter :: invradian=0.017453292
      real*8, parameter :: rearth=6371001.0     ! Radius of earth
      real*8  rlon1,rlat1,rlon2,rlat2           ! Pos. in radians
      real*8  x1,y1,z1,x2,y2,z2                 ! Cartesian position
      real*8  dx,dy,dz,dr,dott                  ! Cartesian distances
c
      rlon1=lon1*invradian             !lon1 in rad
      rlat1=(90.-lat1)*invradian       !90-lat1 in rad 
c
      rlon2=lon2*invradian             !lon2 in rad
      rlat2=(90.-lat2)*invradian       !90-lat2 in rad 
c --- x,y,z of pos 1.
      x1= SIN(rlat1)*COS(rlon1)        
      y1= SIN(rlat1)*SIN(rlon1)
      z1= COS(rlat1) 
c --- x,y,z of pos 2.
      x2= SIN(rlat2)*COS(rlon2)       
      y2= SIN(rlat2)*SIN(rlon2)
      z2= COS(rlat2) 
c --- distances in x, y, z 
      dx=x2-x1                        
      dy=y2-y1
      dz=z2-z1
c --- Final calxulations
      dott=max(-1.,min(x1*x2+y1*y2+z1*z2,1.))
      !dr=SQRT(dx*dx+dy*dy+dz*dz)       !distance pytagaros
      !dr=acos(x1*x2+y1*y2+z1*z2)       ! Arc length
      dr=acos(dott)       ! Arc length
      spherdist=dr*rearth
      end function spherdist
c
c
c --- Locate nearest point with gmsk=true, starting search in 
C --- point (ipiv,jpiv)
      subroutine nearestpoint(glon,glat,nx,ny,lon,lat,ipiv,jpiv,
     &   a1,a2,a3,a4,gmsk,ass)
      use mod_xc
      implicit none

      integer, intent(in)  :: nx,ny
      real,    intent(in)  :: glon(nx,ny),glat(nx,ny)
      real,    intent(in)     :: lon,lat
      integer, intent(inout)  :: ipiv,jpiv
      real,    intent(out)    :: a1,a2,a3,a4
      logical, intent(in)     :: gmsk(nx,ny)
      logical, intent(out)    :: ass

      real, allocatable  :: A(:,:)
      integer isize
      integer ia,ib,ja,jb
      integer iloc(2),i,j

      iloc=-1

      do isize=1,20,2
         ia=max(1,ipiv-isize)
         ib=min(nx,ipiv+isize)
         ja=max(1,jpiv-isize)
         jb=min(ny,jpiv+isize)

         allocate( A(ia:ib,ja:jb) ) 
         
         do j=ja,jb
         do i=ia,ib
            if (gmsk(i,j)) then
               A(i,j)=spherdist(glon(i,j),glat(i,j),lon,lat)
            else
               A(i,j)=1.0E20
            endif
         enddo
         enddo

         if (minval(A) < 1.0E20) then
            iloc=minloc(A)
            deallocate(A)
            exit
         endif
         deallocate(A)
      enddo

      !   print *,'OLD pivots:',ipiv,jpiv
      ipiv=ia+iloc(1)-1
      jpiv=ja+iloc(2)-1
      !   print *,'NEW pivots:',ipiv,jpiv
      if (gmsk(ipiv,jpiv)) then
         a1=1.0
         ass=.true.
      else
         if (mnproc==1)
     &   print *,'nearest_point (WARNING): Could not find (in,jn)'//
     &           'for ipiv, jpiv=',ipiv,jpiv
         a1=0.0
         ass=.false.
      endif
      
      a2=0.0
      a3=0.0
      a4=0.0
      end subroutine nearestpoint

c --- -----------------------------------------------------------------
c --- Replacement of horizontal interpolation routines bilin and bicubic
c --- with a single routine.
c --- 
c --- "interpug" - interpolation routine for "uniformly spaced" (in lon
c --- lat space) and global data sets which are periodic in longitude. 
c --- Interpolates data from data grid (old) to input grid (new) with 
c --- positions specified by newlon, newlat
c --- NB: Also handles Gaussian grids from NCEP. So, not "uniform" in latitude
c ---
c --- Interpolation methods:
c --- itype == 0 : bilinear interpolation. Continous, but discontinuous
c --- 1st derivatives. Always monotonic
c ---
c --- itype == 1 : bicubic interpolation. Continous 0-th and 1st order 
c --- derivatives. 2nd order cross derivative continous at corners. Not 
c --- Monotonic
c ---
c --- Assumptions: 
c ---  1) Data grid is uniform in lon and lat directions. Grid indices
c ---     increase as longitude/latitude increases.
c ---  OR Data grid is uniform in lon direction. Grid indices increase as 
c ---     longitude increases. Latitude uses "Gaussian" points (right now
c ---     any monotonically increasing/decreasing latitude vector will
c ---     work).
c ---  2) Grid is periodic in longitude direction. Should cover the 
c ---     globe.
c --- -----------------------------------------------------------------
      subroutine interpug(old,onx,ony,olonref,olatref,odlon,odlat,
     &                   new,newlon,newlat,itype,gausslat)   
      use mod_xc
      implicit none
c --- Dims of old (data) grid, as well as reference lon/lat and grid increment
      integer, intent(in) :: onx, ony,itype 
      real, intent(in)    :: olonref, olatref    
      real, intent(in)    :: odlon  , odlat
      real, intent(in)    :: old(onx,ony)! old grid
c
c --- Longitude/Latitude of new grid, as well as new (interpolated) values
      real, intent(in), dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &   newlon, newlat
      real, intent(out)   :: new   (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)  
c
c --- Optional argument indicating gaussian grid
      real, intent(in), optional  :: gausslat(ony)
c
      integer i,j,ipos,jpos,ia,ib,ja,jb,ib2,jb2,numerr,j2
      real    aa,bb,lon,minlat,rnumerr, newlon2,dlatg
      real :: f00,f10,f01,f11,fx00,fx10,fx01,fx11,fy00,fy10,fy01,fy11,
     &        fxy00,fxy10,fxy01,fxy11, a1,a2,a3,a4
      real :: rhs(16), coeffs(16), radian
      logical :: gaussgrid
      real, parameter :: thlat=30. ! TODO: tune this parameter
      radian=asin(1.)/90.
c
c --- Security check
      if (itype<0 .and. itype>1) then
         if (mnproc==1) write(lp,'(a)') 'Invalid option for itype'
         call xcstop('(mod_hycom_nersc:interpug)')
      end if
c
c --- Check for gaussian
      gaussgrid=.false.
      if (present(gausslat)) then
         gaussgrid=.true.
      end if
c
c --- Start interpolation
      numerr=0
C$OMP PARALLEL DO PRIVATE(i,j,ipos,jpos,jb,lon,ib, newlon2,
C$OMP& ia,ja, ib2, jb2, f00,f10,f01,f11,fx00,fx10,fx01,fx11,
C$OMP& fy00,fy10,fy01,fy11,fxy00,fxy10,fxy01,fxy11,
C$OMP& a1,a2,a3,a4,rhs,coeffs,j2,dlatg)
C$OMP& SCHEDULE(STATIC,jblk) REDUCTION(+:numerr)
      do j=1-nbdy,jj+nbdy
      do i=1-nbdy,ii+nbdy
c ---    New index in old data. New point is between
c ---    ipos, ib and jpos, jb
cLB         newlon2 = mod(newlon(i,j)+360.d0,360.d0)
         newlon2 = mod(newlon(i,j)-olonref+360.d0,360.d0) + olonref 
         ipos =int((newlon2-olonref)/odlon+1.0d0)  ! [1:onx]
         if (ipos>onx .or. ipos < 1 ) then
            write(lp,'(a,2i5,2f10.2)') 'ipos error ',
     &         ipos,onx,newlon2,olonref+odlon*(onx-1)
            numerr=numerr+1
         endif   
         ia =mod(onx+ipos-2,onx)+1    ! Periodic
         ib =mod(ipos,onx)+1          ! Periodic
         ib2=mod(ib  ,onx)+1          ! Periodic
c
c ---    jpos on uniform grid
         if (.not.gaussgrid) then
            jpos =int((newlat(i,j)-olatref)/odlat+1.0d0)
            if (jpos>ony .or. jpos < 1 ) then
               write(lp,'(a,2i5,2f10.2)') 'jpos error ',
     &            jpos,ony,newlat(i,j),olatref+odlat*(ony-1)
               numerr=numerr+1
            endif   
c ---       TODO: Latitude "Wrap-over" at North Pole ?
            ja   =min(max(1,jpos-1),ony)
            jb   =min(max(1,jpos+1),ony)
            jb2  =min(max(1,jpos+2),ony)
c
c ---    jpos on gaussian grid. NB: There exists formulas for the
c ---    Gaussian locations, and they should be used 
         else
            if (odlat>0.) then ! latitude inrease with increasing index
               jpos=ony
               do j2=1,ony
                  if (gausslat(j2)>newlat(i,j)) then
                     jpos=j2
                     exit
                  end if
               end do
c ---          TODO: Latitude "Wrap-over" at North Pole ?
               jpos =min(max(1,jpos  ),ony)
               ja   =min(max(1,jpos-1),ony)
               jb   =min(max(1,jpos+1),ony)
               jb2  =min(max(1,jpos+2),ony)
            else ! Latitude decreases with increasing index
               jpos=1
               do j2=1,ony
                  if (gausslat(j2)<newlat(i,j)) then
                     jpos=j2
                     exit
                  end if
               end do
c ---          TODO: Latitude "Wrap-over" at North Pole ?
               jpos =min(max(1,jpos  ),ony)
               ja   =min(max(1,jpos+1),ony)
               jb   =min(max(1,jpos-1),ony)
               jb2  =min(max(1,jpos-2),ony)
            end if
         endif   
C        if (j+j0==250) print '(a,5i4)','i',i,ia,ipos,ib,ib2
c
c ---    Grid distance new point -> ipos, jpos. aa,bb in [0,1]
         aa=(newlon2 - olonref-real(ipos-1)*odlon)/odlon
         if (.not. gaussgrid) then
            bb=(newlat(i,j) - olatref-real(jpos-1)*odlat)/odlat
         else
            if (gausslat(jb)==gausslat(jpos)) then
               dlatg=abs(90-gausslat(jpos))
            else
               dlatg=gausslat(jb)-gausslat(jpos)
            end if
            dlatg =abs(dlatg) 
            bb=(newlat(i,j) - gausslat(jpos))/dlatg
         end if
c
c ---    Catch errors - but dont stop until after loop
         if ((aa > 1.0).or.(aa < 0.0)) then
            !write(*,'(3i5,3f10.2)')i,j,ipos,lon,newlon(i,j),lon+odlon
            write(*,'(3i5,3f10.2)')i,j,ipos,lon,newlon2,lon+odlon
            print *,'interpug: invalid aa',aa
            numerr=numerr+1
         endif
         if ((bb > 1.0).or.(bb < 0.0)) then
            if (gaussgrid) then
               write(*,'(3i5,3f10.2)')i,j,jpos,gausslat(jpos),
     &            newlat(i,j),gausslat(jb)
            else
               write(*,'(3i5,3f10.2)')i,j,jpos,lon,newlat(i,j),
     &            olatref+(jpos-1)*odlat
            end if
            print *,'interpug: invalid bb',bb
            numerr=numerr+1
         endif
c
c ---   Set up bilinear weights
        if (itype==0) then
c ---      Bilinear weights
           a1=(1.0-aa)*(1.0-bb)
           a2=aa*(1.0-bb)
           a3=aa*bb
           a4=(1.0-aa)*bb
c ---      New data value
           new(i,j) = a1*old(ipos,jpos)+a2*old(ib  ,jpos)+
     &                a3*old(ib  ,jb  )+a4*old(ipos,jb  ) 
c ---   TODO: most of this can be re-used if (ipos,jpos) hasnt changed
c ---   Set up function and derivatives at ecmwf nodes
c ---   TODO: Change plat in derivatives to data grid lat values
        elseif (itype==1) then
          f00  =old(ipos,jpos)
          f10  =old(ib  ,jpos)
          f01  =old(ipos,jb  )
          f11  =old(ib  ,jb  )
c ---     X derivative with gridspacing 1 - LB: no gridsize needed
          fx00 = 0.5*(old(ib  ,jpos) - old(ia  ,jpos))
          fx10 = 0.5*(old(ib2 ,jpos) - old(ipos,jpos))
          fx01 = 0.5*(old(ib  ,jb  ) - old(ia  ,jb  ))
          fx11 = 0.5*(old(ib2 ,jb  ) - old(ipos,jb  ))
c ---     Y derivative with gridspacing 1 
          fy00 = 0.5*(old(ipos,jb  ) - old(ipos,ja  ))
          fy10 = 0.5*(old(ib  ,jb  ) - old(ib  ,ja  ))
          fy01 = 0.5*(old(ipos,jb2 ) - old(ipos,jpos))
          fy11 = 0.5*(old(ib  ,jb2 ) - old(ib  ,jpos))
c ---     Cross derivative with gridspacing 1  
          fxy00=0.25*( old(ib  ,jb )-old(ib  ,ja  )-
     &                (old(ia  ,jb )-old(ia  ,ja  )))
          fxy10=0.25*( old(ib2 ,jb )-old(ib2 ,ja  )-
     &                (old(ipos,jb )-old(ipos,ja  )))
          fxy01=0.25*( old(ib  ,jb2)-old(ib  ,jpos)-
     &                (old(ia  ,jb2)-old(ia  ,jpos)))
          fxy11=0.25*( old(ib2 ,jb2)-old(ib2 ,jpos)-
     &                (old(ipos,jb2)-old(ipos,jpos)))
c ---     RHS of coeff equation
          rhs=(/f00,f10,f01,f11,fx00,fx10,fx01,fx11,fy00,fy10,fy01,fy11,
     &        fxy00,fxy10,fxy01,fxy11/)
c ---     Solve matrix for cubic coeffs
c ---     TODO: optimize this routine
          coeffs=cubiccoeff(rhs)
c ---     Calculate solution
          new(i,j)=cubicsol(coeffs,aa,bb)
        end if
      end do
      end do
c --- Halt on errors
      rnumerr=numerr ;call xcmaxr(rnumerr); numerr=int(rnumerr)
      if (numerr>0) then
         if (mnproc==1)write(lp,'(a)')'Error(s) occured in interpug..'
         call xcstop('(interpug)')
         stop '(interpug)'
      end if
      end subroutine interpug


c
c --- Routine calculates bilinear weights for interpolation of field
c --- "old" to "new" locations. Also does the interpolation. A Periodic
c --- "old" grid is assumed
c --- TODO: Consider removing this routine
      subroutine bilin_ecmwf2(old,onx,ony,olonref,olatref,odlon,odlat,
     &                       new,newlon,newlat)   
      use mod_xc
      implicit none
      integer, intent(in) :: onx         ! x-dimension of old field
      integer, intent(in) :: ony         ! y-dimension of old field
      real, intent(in)    :: old(onx,ony)! old grid
      real, intent(in)    :: olonref     ! Lon - reference point
      real, intent(in)    :: olatref     ! Lat - reference point
      real, intent(in)    :: odlon       ! Lon - grid spacing in old grid
      real, intent(in)    :: odlat       ! Lat - grid spacing in old grid
c
      real, intent(out)   :: new   (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)  ! New interpolated field
      real, intent(in)    :: newlon(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ! Longitudes for new grid
      real, intent(in)    :: newlat(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ! Latitudes for new grid
c
      integer i,j,ia,ib,ja,jb,ifalse,maxlati,minlati, numerr
      integer ipos    !  index of i-pivot grid point in old grid
      integer jpos    !  index of j-pivot grid point in old grid
      real    aa,bb,a1,a2,a3,a4,maxlat,lon,minlat,rnumerr, newlon2
c
      maxlat=olatref+(ony-1)*odlat
      maxlat=max(olatref,maxlat) ; ! For odlat<0
      minlat=olatref+(ony-1)*odlat
      minlat=min(minlat,olatref)   ! For odlat<0
c
c --- Set "minimumlat" and "maximumlat" index..
      if (olatref<olatref+(ony-1)*odlat) then
         minlati=1
         maxlati=ony
      else
         minlati=ony
         maxlati=1
      end if
c
      numerr=0
c --- Start interpolation
C$OMP PARALLEL DO PRIVATE(i,j,ipos,jpos,jb,lon,ib, newlon2,
c$OMP&                    aa,bb,a1,a2,a3,a4) 
C$OMP& SCHEDULE(STATIC,jblk) REDUCTION(+:numerr)
      do j=1-nbdy,jj+nbdy
      do i=1-nbdy,ii+nbdy
c ---    New index in old data. New point is between
c ---    ipos, ipos+1 and jpos, jpos+1
         newlon2 = mod(newlon(i,j)+360.d0,360.d0)
         !ipos =int((newlon(i,j)-olonref)/odlon+1.0)
         ipos =int((newlon2-olonref)/odlon+1.0d0)
         ib=mod(ipos,onx)+1          ! Periodic
         lon=olonref+(ipos-1)*odlon
         !print *,newlon(i,j),newlon2
c
         jpos =int((newlat(i,j)-olatref)/odlat+1.0d0)
         if (jpos>ony .or. jpos < 1 ) then
            write(lp,'(a,2i5,2f10.2)') 'jpos error ',
     &         jpos,ony,newlat(i,j),olatref+odlat*(ony-1)
            numerr=numerr+1
         endif
         jpos =min(max(1,jpos  ),ony)
         jb   =min(max(1,jpos+1),ony)
         !if (j+j0==250) print *,i,ipos,ib,jpos,newlon2
c
c ---    Grid distance new point -> ipos, jpos
         aa=(newlon2 - olonref-float(ipos-1)*odlon)/odlon
         bb=(newlat(i,j) - olatref-float(jpos-1)*odlat)/odlat
c
c ---    Catch errors - but dont stop until after loop
         if ((aa > 1.0).or.(aa < 0.0)) then
            !write(*,'(3i5,3f10.2)')i,j,ipos,lon,newlon(i,j),lon+odlon
            write(*,'(3i5,3f10.2)')i,j,ipos,lon,newlon2,lon+odlon
            print *,'bilin_ecmwf2: invalid aa',aa
            numerr=numerr+1
         endif
         if ((bb > 1.0).or.(bb < 0.0)) then
            print *,'bilin_ecmwf2: invalid bb',bb
            numerr=numerr+1
         endif

c ---    Bilinear weights
         a1=(1.0-aa)*(1.0-bb)
         a2=aa*(1.0-bb)
         a3=aa*bb
         a4=(1.0-aa)*bb
c ---    New data value
         new(i,j) = a1*old(ipos,jpos)+a2*old(ib  ,jpos)+
     &              a3*old(ib  ,jb  )+a4*old(ipos,jb  ) 
      enddo
      enddo
C$OMP END PARALLEL DO
      rnumerr=numerr
      call xcmaxr_0o(rnumerr)
      numerr=int(rnumerr)

      if (numerr>0) then
         if (mnproc==1) then
            write(lp,'(a)') 'Error(s) occured in bilin_ecmwf2..'
            call flush(lp)
         end if
         call xcstop('(bilin_ecmwf2)')
         stop '(bilin_ecmwf2)'
      end if
      end subroutine bilin_ecmwf2
c
c --- TODO: Consider removing this routine in favor of interpug
      subroutine bicubic(old,onx,ony,olonref,olatref,odlon,odlat,
     &                       new,newlon,newlat)   
      use mod_xc
      implicit none
      integer, intent(in) :: onx         ! x-dimension of old field
      integer, intent(in) :: ony         ! y-dimension of old field
      real, intent(in)    :: old(onx,ony)! old grid
      real, intent(in)    :: olonref     ! Lon - reference point
      real, intent(in)    :: olatref     ! Lat - reference point
      real, intent(in)    :: odlon       ! Lon - grid spacing in old grid
      real, intent(in)    :: odlat       ! Lat - grid spacing in old grid
c
      real, intent(out)   :: new   (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)  ! New interpolated field
      real, intent(in)    :: newlon(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ! Longitudes for new grid
      real, intent(in)    :: newlat(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ! Latitudes for new grid
c
      integer i,j,ia,ib,ja,jb, ib2, jb2
      integer :: maxlati,minlati, numerr
      integer ipos    !  index of i-pivot grid point in old grid
      integer jpos    !  index of j-pivot grid point in old grid
      real    aa,bb,maxlat,lon,minlat,rnumerr, newlon2
      real :: f00,f10,f01,f11,fx00,fx10,fx01,fx11,fy00,fy10,fy01,fy11,
     &        fxy00,fxy10,fxy01,fxy11
      real :: rhs(16), coeffs(16)
c
      maxlat=olatref+(ony-1)*odlat
      maxlat=max(olatref,maxlat) ; ! For odlat<0
      minlat=olatref+(ony-1)*odlat
      minlat=min(minlat,olatref)   ! For odlat<0

      !call xcstop ('(bicubic needs more testing )')     

c
c --- Set "minimumlat" and "maximumlat" index..
      if (olatref<olatref+(ony-1)*odlat) then
         minlati=1
         maxlati=ony
      else
         minlati=ony
         maxlati=1
      end if
c
c
c --- Start interpolation
      numerr=0
C$OMP PARALLEL DO PRIVATE(i,j,ipos,jpos,jb,lon,ib, newlon2,
C$OMP& ia,ja, ib2, jb2, f00,f10,f01,f11,fx00,fx10,fx01,fx11,
C$OMP& fy00,fy10,fy01,fy11,fxy00,fxy10,fxy01,fxy11,
C$OMP& rhs,coeffs)
C$OMP& SCHEDULE(STATIC,jblk) REDUCTION(+:numerr)
      do j=1-nbdy,jj+nbdy
      do i=1-nbdy,ii+nbdy
c ---    New index in old data. New point is between
c ---    ipos, ipos+1 and jpos, jpos+1
         newlon2 = mod(newlon(i,j)+360.d0,360.d0)
         ipos =int((newlon2-olonref)/odlon+1.0d0)
         ia =mod(onx+ipos-2,onx)+1          ! Periodic
         ib =mod(ipos,onx)+1          ! Periodic
         ib2=mod(ib  ,onx)+1          ! Periodic
         lon=olonref+(ipos-1)*odlon
c
         jpos =int((newlat(i,j)-olatref)/odlat+1.0d0)
         if (jpos>ony .or. jpos < 1 ) then
            write(lp,'(a,2i5,2f10.2)') 'jpos error ',
     &         jpos,ony,newlat(i,j),olatref+odlat*(ony-1)
            numerr=numerr+1
         endif
         jpos =min(max(1,jpos  ),ony)
         ja   =min(max(1,jpos-1),ony)
         jb   =min(max(1,jpos+1),ony)
         jb2  =min(max(1,jpos+2),ony)
C         if (j+j0==250) print '(a,5i4)','i',i,ia,ipos,ib,ib2

c
c ---    Grid distance new point -> ipos, jpos
         aa=(newlon2 - olonref-float(ipos-1)*odlon)/odlon
         bb=(newlat(i,j) - olatref-float(jpos-1)*odlat)/odlat
c
c ---    Catch errors - but dont stop until after loop
         if ((aa > 1.0).or.(aa < 0.0)) then
            !write(*,'(3i5,3f10.2)')i,j,ipos,lon,newlon(i,j),lon+odlon
            write(*,'(3i5,3f10.2)')i,j,ipos,lon,newlon2,lon+odlon
            print *,'bilin_ecmwf2: invalid aa',aa
            numerr=numerr+1
         endif
         if ((bb > 1.0).or.(bb < 0.0)) then
            print *,'bilin_ecmwf2: invalid bb',bb
            numerr=numerr+1
         endif
c ---   TODO: most of this can be re-used if (ipos,jpos) hasnt changed
c ---   Set up function and derivatives at ecmwf nodes
        f00  =old(ipos,jpos)
        f10  =old(ib  ,jpos)
        f01  =old(ipos,jb  )
        f11  =old(ib  ,jb  )
c ---   X derivative with gridspacing 1 - TODO: modify with latitude
        fx00 = 0.5*(old(ib  ,jpos) - old(ia  ,jpos))/odlon
        fx10 = 0.5*(old(ib2 ,jpos) - old(ipos,jpos))/odlon
        fx01 = 0.5*(old(ib  ,jb  ) - old(ia  ,jb  ))/odlon
        fx11 = 0.5*(old(ib2 ,jb  ) - old(ipos,jb  ))/odlon
c ---   Y derivative with gridspacing 1 
        fy00 = 0.5*(old(ipos,jb  ) - old(ipos,ja  ))/odlat
        fy10 = 0.5*(old(ib  ,jb  ) - old(ib  ,ja  ))/odlat
        fy01 = 0.5*(old(ipos,jb2 ) - old(ipos,jpos))/odlat
        fy11 = 0.5*(old(ib  ,jb2 ) - old(ib  ,jpos))/odlat
c ---   Cross derivative with gridspacing 1  - TODO: modify with latitude
        fxy00=0.25*(
     &        old(ib ,jb )-old(ib ,ja  )-(old(ia  ,jb )-old(ia  ,ja  ))
     &        )/(odlon*odlat)
        fxy10=0.25*(
     &        old(ib2,jb )-old(ib2,ja  )-(old(ipos,jb )-old(ipos,ja  ))
     &        )/(odlon*odlat)
        fxy01=0.25*(
     &        old(ib ,jb2)-old(ib ,jpos)-(old(ia  ,jb2)-old(ia  ,jpos))
     &        )/(odlon*odlat)
        fxy11=0.25*(
     &        old(ib2,jb2)-old(ib2,jpos)-(old(ipos,jb2)-old(ipos,jpos))
     &        )/(odlon*odlat)
c ---   Testing
c       fx00 =0.
c       fx10 =0.
c       fx01 =0.
c       fx11 =0.
c       fy00 =0.
c       fy10 =0.
c       fy01 =0.
c       fy11 =0.
c       fxy00=0.
c       fxy10=0.
c       fxy01=0.
c       fxy11=0.
c ---   Form rhs
        rhs=(/f00,f10,f01,f11,fx00,fx10,fx01,fx11,fy00,fy10,fy01,fy11,
     &        fxy00,fxy10,fxy01,fxy11/)
c ---   Solve matrix for cubic coeffs
        coeffs=cubiccoeff(rhs)
c ---   Calculate solution
        new(i,j)=cubicsol(coeffs,aa,bb)
      enddo
      enddo
C$OMP END PARALLEL DO
      rnumerr=numerr
      call xcmaxr_0o(rnumerr)
      numerr=int(rnumerr)
      !call xcstop('bicubic')

      if (numerr>0) then
         if (mnproc==1) then
            write(lp,'(a)') 'Error(s) occured in bicubic..'
            call flush(lp)
         end if
         call xcstop('(bicubic)')
         stop '(bicubic)'
      end if
      end subroutine bicubic

      function cubiccoeff(rhs)
      implicit  none
      real, dimension(16), intent(in) :: rhs
      real, dimension(16) :: cubiccoeff
C --- TODO: Matrix is sparse - so there is room for reducing the number
c --- of operations (i.e avoid matmul)
      real, dimension(16*16), parameter ::
     &invcb=(/ 1, 0,-3, 2, 0, 0, 0, 0,-3, 0, 9,-6, 2, 0,-6, 4,
     &         0, 0, 3,-2, 0, 0, 0, 0, 0, 0,-9, 6, 0, 0, 6,-4,
     &         0, 0, 0, 0, 0, 0, 0, 0, 3, 0,-9, 6,-2, 0, 6,-4,
     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9,-6, 0, 0,-6, 4,
     &         0, 1,-2, 1, 0, 0, 0, 0, 0,-3, 6,-3, 0, 2,-4, 2,
     &         0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 3,-3, 0, 0,-2, 2,
     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 3,-6, 3, 0,-2, 4,-2,
     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 2,-2,
     &         0, 0, 0, 0, 1, 0,-3, 2,-2, 0, 6,-4, 1, 0,-3, 2,
     &         0, 0, 0, 0, 0, 0, 3,-2, 0, 0,-6, 4, 0, 0, 3,-2,
     &         0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 3,-2, 1, 0,-3, 2,
     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 2, 0, 0, 3,-2,
     &         0, 0, 0, 0, 0, 1,-2, 1, 0,-2, 4,-2, 0, 1,-2, 1,
     &         0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 2,-2, 0, 0,-1, 1,
     &         0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 2,-1, 0, 1,-2, 1,
     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0,-1, 1/)
      real,parameter,dimension(16,16):: invcb2=reshape(invcb,(/16,16/))
      cubiccoeff=matmul(invcb2,rhs)
      end function

c --- Final calculation of bicubic interpolation
      real function cubicsol(coeffs,aa,bb)
      implicit none
      real, intent(in) :: aa,bb,coeffs(16)
      real :: coeffs2(4,4)
      integer :: i,j
c --- Reshape to c_ij format (see any work on bicubic interpolation)
      coeffs2=reshape(coeffs,(/4,4/))
      cubicsol=0.
      do j=1,4
      do i=1,4
         cubicsol=cubicsol+coeffs2(i,j)*aa**(i-1)*bb**(j-1)
      end do
      end do
      end function cubicsol



c --- Helper routine for creating netcdf diagnostics. 
c --- no record dim assumed. 2D fields only
c --- (hey, thats why its called "draft" - it aint pretty)
      subroutine ncdraft(ncfil,field,vname,iostatus)
      use mod_xc
      use netcdf
      implicit none
      character(len=*)    , intent(in) :: ncfil,iostatus,vname
      real, dimension(:,:), intent(in) :: field
      integer netcdf_error
      netcdf_error=0
      if (mnproc==1) netcdf_error=ncdraft_1(ncfil,field,vname,iostatus)
      call ncerr(netcdf_error)
      end subroutine


c --- Helper routine for creating netcdf diagnostics. 
c --- no record dim assumed. 1 task version
      integer function ncdraft_1(ncfil,field,vname,iostatus)
      use mod_xc
      use netcdf
      implicit none
      character(len=*)    , intent(in) :: ncfil,iostatus,vname
      real, dimension(:,:), intent(in) :: field
      integer :: ncid, dimids(2),sdims(2), dimlen(2), varid,
     &   ndims, nvars,k,k2
      character(len=5) :: ctag
c
      ndims=0
      if (trim(iostatus)=='clobber') then
         ncdraft_1=NF90_create(ncfil,NF90_CLOBBER,ncid)
         if (ncdraft_1/=NF90_NOERR) return
         ndims=0
      else
         ncdraft_1=NF90_open(ncfil,NF90_WRITE,ncid)
         if (ncdraft_1/=NF90_NOERR) return
         ncdraft_1=nf90_inquire(ncid, nDimensions=ndims, 
     &      nVariables=nvars)
         if (ncdraft_1/=NF90_NOERR) return
         ncdraft_1=nf90_redef(ncid)
         if (ncdraft_1/=NF90_NOERR) return
      end if
      sdims(1)=size(field,1)
      sdims(2)=size(field,2)
      !print *,'size(field)',sdims
      dimids=-1
c --- Look for dimensions with same size as input
      do k=1,ndims
         ncdraft_1=nf90_inquire_dimension(ncid, k, len=dimlen(k))
         if (ncdraft_1/=NF90_NOERR) return
         do k2=1,2
            if (dimlen(k)==sdims(k2)) then
               dimids(k2)=k
            end if
         end do
      end do
      !print *,'dimids:',dimids
c --- if any dimids=-1 we have to create them
      do k=1,2
         if (dimids(k)==-1) then
            if (sdims(1)==sdims(2)) then
             write(ctag,'(i5.5)') k
            else
             write(ctag,'(i5.5)') sdims(k)
            end if
            !print *,'ctag:',ctag
            ncdraft_1=nf90_def_dim(ncid, 'dim'//ctag,
     &          sdims(k), dimids(k))
            if (ncdraft_1/=NF90_NOERR) return
         end if
      end do
c --- TODO - check for var name, TODO - transpose
      ncdraft_1=nf90_def_var(ncid,vname,NF90_FLOAT,dimids,varid)
      if (ncdraft_1/=NF90_NOERR) return
      ncdraft_1=nf90_enddef(ncid)
      if (ncdraft_1/=NF90_NOERR) return
      ncdraft_1=nf90_put_var(ncid,varid,field)
      if (ncdraft_1/=NF90_NOERR) return
      ncdraft_1=nf90_close(ncid)
      if (ncdraft_1/=NF90_NOERR) return
      end function




c --- Routines checks error value of netcdf operation. Exits
c --- if an error occured
      subroutine ncerr(errcode)
      use mod_xc 
      use netcdf
      implicit none
      integer, intent(in) :: errcode
      integer :: errcode2
      real a(1)
      a(1)=errcode
      call xcastr(a,1)
      errcode2=int(a(1))
      if (errcode2/=NF90_NOERR) then
         if (mnproc==1) then
            write(lp,'(a)') NF90_STRERROR(errcode2)
            print *,errcode2
         end if
         stop '(ncerr)'
         call xcstop('(ncerr)')
      end if
      end subroutine
      

      end module
