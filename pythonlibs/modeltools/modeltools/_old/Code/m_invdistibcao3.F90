module m_invdistibcao3
contains
subroutine invdistibcao3(deep,modlon,modlat,nocean,nland,nx,ny,&
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

   real aa,bb,a1,a2,a3,a4,s1,s2,s3,s4
   integer i,j,ipiv,jpiv,ipiv1,jpiv1,ia,ja
   real dxi,dyi
   logical lpiv
   integer minlocd(2)
   real d(1:3,1:3),alpha
   real a(1:3,1:3),suma
   integer ifile
   character(len=3) tag3

   ifile=0

!$OMP PARALLEL DO PRIVATE(i,j,ipiv,jpiv,lpiv,ja,ia,d,minlocd,suma,a) SHARED(nxd,nyd) SCHEDULE(STATIC,1)
   do j=1,ny
   ipiv=nxd/2
   jpiv=nyd/2
   do i=1,nx
      if ((minlon  < modlon(i,j) .and. modlon(i,j) < maxlon + dx).and.&
          (minlat  < modlat(i,j) .and. modlat(i,j) < maxlat + dy)) then

         call pivotsearch(modlon(i,j),modlat(i,j),lon,lat,nxd,nyd,ipiv,jpiv,lpiv)

         if (lpiv .and. 1<ipiv .and. ipiv<nxd-2 .and.  1<jpiv .and. jpiv<nyd-2) then

            do ja=1,3
            do ia=1,3
               d(ia,ja)=spherdist(modlon(i,j),modlat(i,j),&
                                  lon(ipiv+ia-2,jpiv+ja-2),lat(ipiv+ia-2,jpiv+ja-2))
            enddo
            enddo

            if (minval(d)/maxval(d) < 0.00001) then
               minlocd=minloc(d)
               deep(i,j)=obs(ipiv-2+minlocd(1),jpiv-2+minlocd(2))
            else

               suma=0.0
               do ja=1,3
               do ia=1,3
                  a(ia,ja)=1.0/d(ia,ja)**3.5
                  suma=suma+a(ia,ja)
               enddo
               enddo

               deep(i,j)=0.0
               do ja=1,3
               do ia=1,3
                  deep(i,j)=deep(i,j)+a(ia,ja)*obs(ipiv+ia-2 ,jpiv+ja-2 )
               enddo
               enddo
               deep(i,j)=deep(i,j)/suma
            endif

         endif
      endif

   enddo
   enddo
end subroutine invdistibcao3
end module m_invdistibcao3
