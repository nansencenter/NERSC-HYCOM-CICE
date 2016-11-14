module mod_conman
! IBCAO variables
integer, parameter :: conman_nx=720
integer, parameter :: conman_ny=578

type CONMAN_data
   real deep  (conman_nx,conman_ny)
   real lat   (conman_nx,conman_ny)
   real lon   (conman_nx,conman_ny)
   real weight(conman_nx,conman_ny)
   real minlon,maxlon,minlat,maxlat
end type CONMAN_data
type (CONMAN_data) conman


contains 
subroutine read_felt()
IMPLICIT NONE

character(len=200) :: path0
CHARACTER(LEN=*),parameter  :: ifile1=  'CONMAN4km-3.seq_sgi'

! Input file with model data
CHARACTER(LEN=99)  :: ofile        ! Output file

INTEGER, PARAMETER :: im=720, jm=578 ! Grid dimensions (alternativ, read std. inp. and allocate arrays)
REAL               :: field(im,jm) ! Field
INTEGER            :: mode         ! Choose how to read felt-file
REAL               :: gparam(6)    ! Grid parameters
INTEGER            :: ident(20)    ! Identity string for the sequential data
INTEGER            :: igtype       ! Grid type
INTEGER            :: ierr         ! Error indicator

real :: lon,lat,dx,dy,dz,rp,pi,om,an,fi
integer :: i,j
character(len=230) :: ifile
logical :: ex


! Retrieve etopo5 path
call getenv('CONMAN_PATH',path0)
if (trim(path0)=='') then
   print *,'No path set for CONMAN. Set the environment variable CONMAN_PATH'
   print *,'To point to the location of the CONMAN bathymetry'
   call exit(1)
end if
path0=trim(path0)//'/'


ifile=trim(path0)//trim(ifile1)


inquire(file=trim(ifile),exist=ex)
if (.not.ex) then
   print *,'Can not find CONMAN file '//trim(ifile)
   print *,'path0=',trim(path0)
   stop '(mod_conman)'
end if

OPEN(UNIT=11,FILE=trim(ifile),ACCESS='SEQUENTIAL',STATUS='UNKNOWN',FORM='UNFORMATTED')
mode = 0   ! Read fields sequentially
CALL getflt(mode,11,0,ident,im*jm,field,igtype,gparam,0,ierr)
IF (ierr /= 0) THEN
  PRINT *,"Error with reading file ",ifile
  STOP
ENDIF
CLOSE(11)
print*, 'CONMAN Ident is  ', ident
!print*, 'CONMAN Topo Test ', field(360,:)
print*, 'CONMAN gparam is  ', gparam(:)
conman%deep = field


conman%minlat= 400
conman%maxlat=-400
conman%minlon= 400
conman%maxlon=-400

! Convert felt polar stereographic coordinates to lon/lat
pi = 4.*ATAN(1.)
om = 180./pi
AN = (79.*150)/(ABS(ident(17)*0.1))
print *,'an1:',an
AN = gparam(3)
print *,'an2:',an
!fi=0.
fi=gparam(4)
print *,'fi:',fi
DO j=1,jm
DO i=1,im
   dx = REAL(i) - gparam(1)
   dy = gparam(2) - REAL(j)
   rp = (dx*dx+dy*dy)**0.5
   lat = 90.-om*2.*ATAN(rp/an)
   lon = 0.
   IF (rp  > 1.E-10) lon = fi + om*ATAN2(dx,dy)
   IF (lon <= -180.) lon = 360. + lon
   IF (lon >   180.) lon = lon - 360.
   conman%lat(i,j)=lat
   conman%lon(i,j)=lon

   ! Linear weight
   conman%weight(i,j)=min(                                 &
                       min(float(          i)/20.0, 1.0),  &
                       min(float(          j)/20.0, 1.0),  &
                       min(float(conman_nx-i)/20.0, 1.0),  &
                       min(float(conman_ny-j)/20.0, 1.0) )

   conman%maxlat=max(conman%maxlat,lat)
   conman%minlat=min(conman%minlat,lat)
   conman%maxlon=max(conman%maxlon,lon)
   conman%minlon=min(conman%minlon,lon)

   conman%deep(i,j)=max(0.,conman%deep(i,j))
ENDDO
ENDDO





END subroutine read_felt
end module mod_conman

