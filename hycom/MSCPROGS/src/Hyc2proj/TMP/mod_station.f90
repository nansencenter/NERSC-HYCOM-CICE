










module mod_station


   !KAL  Format of the infile: Example 1  with 3 sections - 
   !KAL  station_group name interpreted as section name

   !#Station_group_name_1
   !lon lat
   !lon lat
   !lon lat
   !lon lat
   !#Station_group_name_2
   !lon lat
   !lon lat
   !lon lat
   !lon lat
   !#Station_group_name_3
   !lon lat
   !lon lat
   !lon lat
   !lon lat


   !KAL  Format of the infile: Example 1  with 3 moorings - 
   !KAL  station_group name interpreted as mooring name

   !#Station_group_name_1
   !lon lat
   !#Station_group_name_2
   !lon lat
   !#Station_group_name_3
   !lon lat


   character(len=*), parameter  :: infile_stations   = 'stations.in'

   type station
      real lon
      real lat
      integer ipiv
      integer jpiv
      real a(4)
      real depth ! water depth
      LOGICAL :: FLAG
   end type station

   type(station), save,  allocatable, dimension(:,:) :: stations
   character(len=40) , save, dimension(:), allocatable :: groupname
   integer , save :: ngroup, maxstat
   integer , save, dimension(:), allocatable :: stations_per_group


   real, parameter, private :: undef=-1e14
contains

subroutine get_stations()
   use mod_grid, only : plon, plat
   use mod_xc  , only : idm, jdm
   use mod_confmap
   implicit none
   logical :: ex
   integer :: ios,counter,prevsec
   integer :: nstat,ipiv,jpiv
   character(len=80) :: c80
   real :: lon,lat,lon_n,lat_n,a(4)

   ! Get stations to process
   inquire(exist=ex,file=infile_stations)
   if (.not. ex) then
      print '(a)','You must specify stations to use in'
      print '(a)','the file called  '//infile_stations
      stop '(get_stations)'
   endif

   ! Start reading to get a) number of groups b) max number of stations per group
   open(10,file=infile_stations,form='formatted')
   ios=0
   ngroup=0
   nstat=0
   maxstat=0
   read(10,'(a)',iostat=ios) c80
   do while (ios==0)
      !print *,'|'//trim(c80)//'|'
      if (c80(1:1)=='#') then
         if (ngroup/=0) then
            maxstat=max(maxstat,nstat)
         end if
         nstat=0
         ngroup=ngroup+1
      else
         nstat=nstat+1
      end if
      read(10,'(a80)',iostat=ios) c80
   end do
   print *,ngroup, nstat
   maxstat=max(maxstat,nstat)
   close(10)

   if (maxstat==0 .or. ngroup==0)  then
      print *,'No  stations or no groups'
      print *,'Groups      :',ngroup
      print *,'Max stations:',maxstat
      stop '(get_stations)'
   end if

   call initconfmap(idm,jdm)

   ! Allocate 
   allocate(stations_per_group(ngroup))
   allocate(groupname         (ngroup))
   allocate(stations          (ngroup,maxstat))

   ! Start reading (again)
   open(10,file=infile_stations,form='formatted')
   !rewind(10)
   ios=0
   ngroup=0
   nstat=0
   read(10,'(a80)',iostat=ios) c80
   !print *,ios
   do while (ios==0)
      if (c80(1:1)=='#') then
         nstat=0
         ngroup=ngroup+1
         groupname(ngroup)=c80(2:len_trim(c80))
      else 
         read(c80,*) lon,lat
         nstat=nstat+1
         !print *,ngroup,nstat,lon,lat
         stations_per_group(ngroup)=nstat
         stations(ngroup,nstat)%lon=lon
         stations(ngroup,nstat)%lat=lat

         ! corresponding pivot points
         call oldtonew(lat,lon,lat_n,lon_n)
         call pivotp(lon_n,lat_n,ipiv,jpiv)

         stations(ngroup,nstat)%ipiv=ipiv
         stations(ngroup,nstat)%jpiv=jpiv
         if (  stations(ngroup,nstat)%ipiv >= idm .or. stations(ngroup,nstat)%ipiv < 1 .or. &
               stations(ngroup,nstat)%jpiv >= jdm .or. stations(ngroup,nstat)%jpiv < 1 ) then
            stations(ngroup,nstat)%flag=.false.
            a(1)=1.; a(2:4)=0.
         else
            stations(ngroup,nstat)%flag=.true.
            ! Bilinear interpolation factors
            call bilincoeff(plon,plat,idm,jdm,lon,lat,ipiv,jpiv,&
                             a(1),a(2),a(3),a(4))
         end if
         stations(ngroup,nstat)%a=a
      end if
      read(10,'(a80)',iostat=ios) c80
   end do
   close(10)
end subroutine get_stations

real function fld2d_to_station(fld2d,station2)
use mod_xc , only : idm, jdm
implicit none
real, intent(in) , dimension(idm,jdm) :: fld2d
type(station), intent(in) :: station2

integer :: i,j

! Simple closest point for now
i=station2%ipiv
j=station2%jpiv
if (station2%flag) then
   fld2d_to_station = fld2d(i,j)
else
   fld2d_to_station = undef
end if
end function


end module mod_station

