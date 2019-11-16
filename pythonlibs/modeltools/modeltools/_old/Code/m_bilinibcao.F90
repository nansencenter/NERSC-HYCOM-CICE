module m_bilinibcao
contains
subroutine bilinibcao(deep,modlon,modlat,nocean,nland,nx,ny,&
                 obs,lat,lon,nxd,nyd,minlon,minlat,maxlon,maxlat,dx,dy)
! inverse distance interpolation from obscured IBCAO grid
   use m_pivotsearch
   use m_spherdist
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

   real aa,bb,s1,s2,s3,s4
   integer i,j,ipiv,jpiv,ipiv1,jpiv1,ia,ja
   real dxi,dyi
   logical lpiv
   integer minlocd(2)
   real d1,d2,d3,d4
   real s12,s13,s14
   real s21,s23,s24
   real s31,s32,s34
   real s41,s42,s43
   real p1,p2,p3,p4
   real a1,a2,a3,a4
   real :: en=1.0

   ipiv=nxd/2
   jpiv=nyd/2
   do i=1,nx
   do j=1,ny
      if ((minlon  < modlon(i,j) .and. modlon(i,j) < maxlon + dx).and.&
          (minlat  < modlat(i,j) .and. modlat(i,j) < maxlat + dy)) then

         call pivotsearch(modlon(i,j),modlat(i,j),lon,lat,nxd,nyd,ipiv,jpiv,lpiv)

         if (lpiv .and. 1<ipiv .and. ipiv<nxd-2 .and.  1<jpiv .and. jpiv<nyd-2) then

            p1=obs(ipiv ,jpiv)
            p2=obs(ipiv+1 ,jpiv)
            p3=obs(ipiv+1 ,jpiv+1)
            p4=obs(ipiv ,jpiv+1)

            d1=spherdist(modlon(i,j),modlat(i,j),lon(ipiv,jpiv),lat(ipiv,jpiv))
            d2=spherdist(modlon(i,j),modlat(i,j),lon(ipiv+1,jpiv),lat(ipiv+1,jpiv))
            d3=spherdist(modlon(i,j),modlat(i,j),lon(ipiv+1,jpiv+1),lat(ipiv+1,jpiv+1))
            d4=spherdist(modlon(i,j),modlat(i,j),lon(ipiv,jpiv+1),lat(ipiv,jpiv+1))

            s12=spherdist(lon(ipiv,jpiv),lat(ipiv,jpiv),lon(ipiv+1,jpiv),lat(ipiv+1,jpiv))
            s13=spherdist(lon(ipiv,jpiv),lat(ipiv,jpiv),lon(ipiv+1,jpiv+1),lat(ipiv+1,jpiv+1))
            s14=spherdist(lon(ipiv,jpiv),lat(ipiv,jpiv),lon(ipiv,jpiv+1),lat(ipiv,jpiv+1))

            s21=s12
            s23=spherdist(lon(ipiv+1,jpiv),lat(ipiv+1,jpiv),lon(ipiv+1,jpiv+1),lat(ipiv+1,jpiv+1))
            s24=spherdist(lon(ipiv+1,jpiv),lat(ipiv+1,jpiv),lon(ipiv,jpiv+1),lat(ipiv,jpiv+1))

            s31=s13
            s32=s23
            s34=spherdist(lon(ipiv+1,jpiv+1),lat(ipiv+1,jpiv+1),lon(ipiv,jpiv+1),lat(ipiv,jpiv+1))

            s41=s14
            s42=s24
            s43=s34


            a1=(en*d2*d3*d4)/(s12*s13*s14)
            a2=(d1*en*d3*d4)/(s21*s23*s24)
            a3=(d1*d2*en*d4)/(s31*s32*s34)
            a4=(d1*d2*d3*en)/(s41*s42*s43)

            deep(i,j)=a1*p1+a2*p2+a3*p3+a4*p4

         endif
      endif

   enddo
   enddo
end subroutine bilinibcao
end module m_bilinibcao
