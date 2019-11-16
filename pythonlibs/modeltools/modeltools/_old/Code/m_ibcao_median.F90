module m_ibcao_median
contains
subroutine ibcao_median(deep,modlon,modlat,nocean,nland,nx,ny,&
                 obs,lat,lon,nxd,nyd,minlon,minlat,maxlon,maxlat,dx,dy)
! Interpolate the IBCAO data set to the model grid
!Francois Counillon 28/05/2013
   use m_pivotsearch
   use m_spherdist
   use m_quick_sort
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
   real,    intent(in)    :: dx                ! data dlon (1/12 for ibcao)
   real,    intent(in)    :: dy                ! data dlon (1/12 for ibcao)

   real aa,bb,a1,a2,a3,a4,s1,s2,s3,s4
   integer i,j,ipiv,jpiv,ipiv1,jpiv1,ia,ja
   real dxi,dyi
   real, allocatable ::  tmp(:)
   logical lpiv
   integer npoint,cnt
   real meandx
   character(len=3) tag3
   !compute the grid size in the middle of the grid and take the median among all
   !the wet point.
   meandx=spherdist(modlon(nx/2,ny/2),modlat(nx/2,ny/2),modlon(nx/2,ny/2)+1,modlat(nx/2,ny/2))/1000
   !taking the median of point within [-npoint:npoint]
   npoint=floor(meandx/dx)/2
   print * , 'Fanf test npoint',npoint,meandx,dx
   allocate(tmp((2*npoint+1)*(2*npoint+1)))
   do j=1,ny
   ipiv=nxd/2
   jpiv=nyd/2
   do i=1,nx
      if ((minlon  < modlon(i,j) .and. modlon(i,j) < maxlon + dx).and.&
          (minlat  < modlat(i,j) .and. modlat(i,j) < maxlat + dy)) then
         !gradiant search towards the minimum dist. Efficient as initial point
         !are sorted
         call pivotsearch(modlon(i,j),modlat(i,j),lon,lat,nxd,nyd,ipiv,jpiv,lpiv)
         if (lpiv .and. 1<=ipiv .and. ipiv<=nxd .and.  1<=jpiv .and. jpiv<=nyd) then
               cnt=0
               tmp(:)=999999.
               do ja=max(1,jpiv-npoint),min(jpiv+npoint,nyd)
               do ia=max(1,ipiv-npoint),min(ipiv+npoint,nxd)
                  if (obs(ia,ja)>0) then
                    cnt=cnt+1
                    tmp(cnt)=obs(ipiv,jpiv)
                  endif
               enddo
               enddo
               if (cnt>0)then
                  call QsortC(tmp)
                  if ( mod(cnt, 2) == 0 ) then
                   deep(i,j) = (tmp(cnt/2+1) + tmp(cnt/2))/2.0
                  else
                  deep(i,j) = tmp(cnt/2+1)
                  end if
               else
                  deep(i,j)=0.
               endif
         endif
      endif

   enddo
   print *,j
   enddo
   deallocate(tmp)
   print *,'Finished interpolation'
end subroutine ibcao_median
end module m_ibcao_median
