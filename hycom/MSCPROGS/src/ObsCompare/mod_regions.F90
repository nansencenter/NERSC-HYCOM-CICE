module mod_regions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Format of region file regiondefs.in
! 2                                # Number of regions
! 4 TEST1                          # Number of points and name of 1st region
! -30 40                           # lon/lat pair - 1st point
!  30 40                           # lon/lat pair - 2nd point
!  30 60                           # lon/lat pair - 3rd point
! -30 60                           # lon/lat pair - 4th point
! 4 TEST2                          # Number of points ansd name of 2nd region
! -30 20                           # etc etc ...
!  30 20
!  30 60
! -30 60
!
!
! Methods:
! read_regions - initializes regions from regiondefs.in - CALLED FIRST!
! getlons      - retrieves longitudes of a region
! getlats      - retrieves latitudes of a region
! getnregions  - retrieves number of regions
! getnpoints   - retrieves number of points for a region
! getname      - retrieves name of a region
! getmask      - creates   mask defining a region
! getmask has module procedures (not called directly):
!     -> getmask_llarray - creates mask when lon/lat are arrays
!     -> getmask_llvector- creates mask when lon/lat are vectors
! dumprms_ascii- dumps values (typically RMS, bias etc) to ascii file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
character(len=*), parameter  :: infile_RMS      = 'regiondefs.in'
integer, parameter :: maxpoints =10
type region
   character(len=20) :: name
   real              :: lons(maxpoints)
   real              :: lats(maxpoints)
   integer           :: npoints
end type region
type(region), allocatable  :: regions(:)
integer, save :: numregion

! NB - private, dont  access directly, but use methods listed above
private :: numregion, regions, region, infile_RMS, maxpoints 

! Interface to handle two typical situations
interface getmask
module procedure     &
   getmask_llarray , &
   getmask_llvec   , &
   getmask_llscalar
end interface
   
contains

! Sets up private variables for this module
subroutine read_regions()
   implicit none

   character(len=20) :: tmpchar

   logical :: ex
   integer :: ios, region_counter,numpoints, point_counter
   real :: v1,v2


   print *

   ! This specifies regions for computing the RMS over
   inquire(exist=ex,file=infile_RMS)
   if (.not. ex) then
      print '(a)','You must specify regions to use in'
      print '(a)','the file called  '//infile_RMS
      stop '(read_regions)'
   endif

   ! Start reading
   ios=0
   region_counter=1
   open(10,file=infile_RMS,form='formatted')
   read(10,*,iostat=ios) numregion
   print '(i3,a)',numregion,' Regions:'
   allocate(regions(numregion))
   do while (region_counter<=numregion .and. ios==0 )

      read (10,*,iostat=ios) numpoints, tmpchar
      if (ios==0) then
         regions(region_counter)%name=trim(tmpchar)

         if (numpoints>maxpoints) then
            print *,'To many points for region number ',region_counter
            stop '(get_regions)'
         end if

         if (numpoints<3) then
            print *,'To few points for region number ',region_counter
            stop '(get_regions)'
         end if

         regions(region_counter)%npoints=numpoints


         write(*,*) 'Number of points in region no ',region_counter,': ',numpoints
         point_counter=1
         do while (point_counter<=numpoints .and. ios==0)

            read (10,*,iostat=ios) v1,v2
            if (ios==0) then
               regions(region_counter)%lons(point_counter)=v1
               regions(region_counter)%lats(point_counter)=v2

               print '(a,2f10.2)', 'lon lat:', &
                        regions(region_counter)%lons(point_counter), &
                        regions(region_counter)%lats(point_counter)
               point_counter=point_counter+1
            else
               print *,'An error occured when reading reguion info'
               stop '(get_regions)'
            end if
         end do

      end if

      region_counter=region_counter+1
   end do

   if (ios/=0) then 
      print *,'An error occured when reading reguion info'
      stop '(get_regions)'
   end if
            

end subroutine read_regions

! Retrieves longitudes of a region
function getlons(iregion)
   implicit none
   integer, intent(in) :: iregion
   real, allocatable, dimension(:) :: getlons
   if (allocated(regions)) then
      if (iregion>0 .and. iregion <=numregion) then
         allocate(getlons(regions(iregion)%npoints))
         getlons=regions(iregion)%lons(1:regions(iregion)%npoints)
      end if
   end if
end function


! Retrieves latitudes of a region
function getlats(iregion)
   implicit none
   integer, intent(in) :: iregion
   real, allocatable, dimension(:) :: getlats
   if (allocated(regions)) then
      if (iregion>0 .and. iregion <=numregion) then
         allocate(getlats(regions(iregion)%npoints))
         getlats=regions(iregion)%lats(1:regions(iregion)%npoints)
      end if
   end if
end function


! Retrieves number of points defined for a region
function getnpoints(iregion)
   implicit none
   integer, intent(in) :: iregion
   integer :: getnpoints
   getnpoints=-1
   if (allocated(regions)) then
      getnpoints=regions(iregion)%npoints
   end if
end function

! Retrieves number of regions
function getnregions()
   implicit none
   integer :: getnregions
   getnregions=-1
   if (allocated(regions)) then
      getnregions=numregion
   end if
end function

! Retrieves region name
function getname(iregion)
   implicit none
   integer, intent(in) :: iregion
   character(len=20) :: getname
   getname=''
   if (allocated(regions)) then
      getname=regions(iregion)%name
   end if
end function

! Retrieves mask when lon/lat are arrays
subroutine getmask_llarray(iregion,mask,modlon,modlat,idm,jdm)
   use mod_sphere_tools
   implicit none
   integer, intent(in)  :: iregion,idm,jdm
   logical, intent(out) :: mask(idm,jdm)
   real   , intent(in ) :: modlon(idm,jdm)
   real   , intent(in ) :: modlat(idm,jdm)
   integer :: i,j,npoints
   mask=.false.
   if (allocated(regions)) then
      npoints=regions(iregion)%npoints
      do j=1,jdm
      do i=1,idm
         mask(i,j)=inbox(regions(iregion)%lons(1:npoints), &
                         regions(iregion)%lats(1:npoints), &
                         npoints,                          &
                         modlon(i,j),modlat(i,j))
      end do
      end do
   end if
end subroutine


! Retrieves mask when lon/lat are scalars
subroutine getmask_llscalar(iregion,mask,modlon,modlat)
   use mod_sphere_tools
   implicit none
   integer, intent(in)  :: iregion
   logical, intent(out) :: mask
   real   , intent(in ) :: modlon
   real   , intent(in ) :: modlat
   integer :: i,j,npoints
   mask=.false.
   npoints=regions(iregion)%npoints
   mask=inbox(regions(iregion)%lons(1:npoints), &
              regions(iregion)%lats(1:npoints), &
              npoints,                          &
              modlon,modlat)
end subroutine

! Retrieves mask when lon/lat are vectors
subroutine getmask_llvec(iregion,mask,modlon,modlat,idm,jdm)
   use mod_sphere_tools
   implicit none
   integer, intent(in)  :: iregion,idm,jdm
   logical, intent(out) :: mask(idm,jdm)
   real   , intent(in ) :: modlon(idm)
   real   , intent(in ) :: modlat(jdm)
   integer :: i,j,npoints
   mask=.false.
   if (allocated(regions)) then
      npoints=regions(iregion)%npoints
      do j=1,jdm
      do i=1,idm
         mask(i,j)=inbox(regions(iregion)%lons(1:npoints), &
                         regions(iregion)%lats(1:npoints), &
                         npoints,                          &
                         modlon(i),modlat(j))
      end do
      end do
   end if
end subroutine


subroutine dumprms_ascii(iregion,yfrac,meanmod,meanobs,bias,rms,npoint,prefix)
real,    intent(in) :: yfrac,meanmod,meanobs,bias,rms
integer, intent(in) :: iregion, npoint
character(len= *), intent(in) :: prefix
character(len=20) :: regname
regname=getname(iregion)
open(11,file=trim(prefix)//trim(regname)//'.asc', position='append', status='unknown')
write(11,'(5f14.5,i8)') yfrac, meanmod, meanobs,bias,rms,npoint
close(11)
end subroutine


end module mod_regions
