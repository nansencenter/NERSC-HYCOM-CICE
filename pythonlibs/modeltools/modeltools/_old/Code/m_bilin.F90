module m_bilin
contains
subroutine bilin(deep,modlon,modlat,nocean,nland,nx,ny,&
                 obs,lat,lon,nxd,nyd,minlon,minlat,maxlon,maxlat,dx,dy,periodic)
! bilinear interpolation from regular lon-lat grid
   implicit none
   integer, intent(in)    :: nx                ! model grid i-dimension
   integer, intent(in)    :: ny                ! model grid j-dimension
   real,    intent(inout) :: deep(nx,ny)       ! Computed depth
   real,    intent(in)    :: modlon(nx,ny)     ! model longitudes
   real,    intent(in)    :: modlat(nx,ny)     ! model latitudes
   integer, intent(in)    :: nocean(nx,ny)     ! Computed depth
   integer, intent(in)    :: nland(nx,ny)      ! Computed depth
   integer, intent(in)    :: nxd               ! data grid i-dimension
   integer, intent(in)    :: nyd               ! data grid j-dimension
   real,    intent(in)    :: obs(nxd,nyd)      ! Observed depth
   real,    intent(in)    :: lat(nxd,nyd)      ! data longitudes
   real,    intent(in)    :: lon(nxd,nyd)      ! data latitudes
   real,    intent(in)    :: minlon            ! data minimum longitude
   real,    intent(in)    :: minlat            ! data minimum latitudes
   real,    intent(in)    :: maxlon            ! data maximum longitude
   real,    intent(in)    :: maxlat            ! data maximum latitudes
   real,    intent(in)    :: dx                ! data dlon (1/12 for etopo5)
   real,    intent(in)    :: dy                ! data dlat (1/12 for etopo5)
   logical, intent(in)    :: periodic          ! Allows for periodic global data sets

   real s1,s2,s3,s4
   integer i,j,ipiv,jpiv,ipiv1,jpiv1
   real*8 dxi,dyi,tmplon,tmplat,aa,bb,a1,a2,a3,a4

   dxi=1.0/dx
   dyi=1.0/dy
   print *,'bilin:',minlon,minval(modlon)

   do i=1,nx
   do j=1,ny

      tmplon=modlon(i,j)
      tmplat=modlat(i,j)

      ! temporary stuff - mostly a issue  on very fine grids
      if (tmplon<=minlon) tmplon=tmplon+360.

      if ((minlon  < tmplon .and. tmplon <= maxlon + dx).and.&
          (minlat  < tmplat .and. tmplat <= maxlat + dy)) then

         ipiv=int((tmplon-minlon)*dxi+1.0)
         jpiv=int((tmplat-minlat)*dyi+1.0)

         if (1 > ipiv .or. ipiv > nxd .or. 1 > jpiv .or. jpiv > nyd ) then
            print *,'bilin: ipiv or jpiv problem: ',ipiv,jpiv,nxd,nyd
            print *,'lon/lat',tmplon,tmplat,minlon,minlat
            print *,dx,dy
            stop
         endif

         if (periodic) then
            ipiv1=mod(ipiv,nxd)+1
         else
            ipiv1=min(ipiv+1,nxd)
         endif
         jpiv1=min(jpiv+1,nyd)

         aa=(tmplon-lon(ipiv,jpiv))*dxi
         bb=(tmplat-lat(ipiv,jpiv))*dyi

         aa=min(1.,max(0.,aa))
         !if (aa<0.or.aa>1) print *,i,j,aa



         a1=(1.0-aa)*(1.0-bb)
         a2=aa*(1.0-bb)
         a3=aa*bb
         a4=(1.0-aa)*bb

         deep(i,j) = a1*obs(ipiv,jpiv)  +a2*obs(ipiv1,jpiv)&
                    +a3*obs(ipiv1,jpiv1)+a4*obs(ipiv,jpiv1)

         if (i==276 .and. j==576) then
            print *,'test:',i,j,aa,deep(i,j)
         end if
      else
         if (i==276 .and. j==576) then
            print *,'test2:',i,j,minlon,tmplon,maxlon+dx
         end if
      endif

   enddo
   enddo
end subroutine bilin
end module m_bilin
