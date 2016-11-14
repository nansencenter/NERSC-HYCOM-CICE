module m_pivotsearch
contains
subroutine pivotsearch(lon,lat,gridlon,gridlat,nxd,nyd,ipiv,jpiv,lpiv)
   use m_spherdist
   implicit none
   integer, intent(in)    :: nxd               ! data grid i-dimension
   integer, intent(in)    :: nyd               ! data grid j-dimension
   real,    intent(in)    :: gridlon(nxd,nyd)  ! data longitudes
   real,    intent(in)    :: gridlat(nxd,nyd)  ! data latitudes
   real,    intent(in)    :: lat               ! longitude location
   real,    intent(in)    :: lon               ! latitude location
   integer, intent(inout) :: ipiv              ! i-pivot point
   integer, intent(inout) :: jpiv              ! j-pivot point
   logical, intent(out)   :: lpiv              ! true if pivot point is found

   real dist(-1:1,-1:1)
   real corner(2,2)
   integer pivotp(2)
   integer mloc(2)
   integer mcor(2)
   integer :: pivcor(2)=2
   integer i,j,k

   pivotp(1)=ipiv
   pivotp(2)=jpiv

   do k=1,10000
      do j=-1,1
      do i=-1,1
         dist(i,j)=spherdist(lon,lat,gridlon(pivotp(1)+i,pivotp(2)+j),&
                                     gridlat(pivotp(1)+i,pivotp(2)+j))
      enddo
      enddo
      mloc=minloc(dist)
      pivotp=pivotp+mloc-pivcor

      if (pivotp(1) < 3 .or. pivotp(1) > nxd-2 .or. &
          pivotp(2) < 3 .or. pivotp(2) > nyd-2) then
          lpiv=.false.
          exit
      endif

      if (mloc(1)==2 .and. mloc(2)==2) then
         corner(1,1)=dist(-1,-1)*dist(-1,0)*dist(0,-1)
         corner(1,2)=dist(-1, 1)*dist(-1,0)*dist(0, 1)
         corner(2,2)=dist( 1, 1)*dist( 1,0)*dist(0, 1)
         corner(2,1)=dist( 1,-1)*dist( 1,0)*dist(0,-1)
!         print '(a,i7,a,2i7)','Iteration=',k,' pivotp=',pivotp
!         print '(3g14.6)',dist
         mcor=minloc(corner)
         where (mcor == 2) mcor=0
         pivotp=pivotp-mcor
         ipiv=pivotp(1)
         jpiv=pivotp(2)
         lpiv=.true.
!         print *,'final pivot point is: ',pivotp
         exit
      endif
   enddo

end subroutine pivotsearch
end module m_pivotsearch
