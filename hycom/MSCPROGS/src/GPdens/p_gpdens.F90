! --- ---------------------------------------------------
! --- Simple program for reading generic gp series
! --- Files needed to run this program are:
! ---  1) a gp-unformatted file (created by HYCOM)
! ---  2) a gp header file      (created by HYCOM)
! ---
! --- Header file contains info on what data can be found 
! --- in one gp-record. This makes this routine 
! --- independent of model size erc...
! ---
! --- Currently extracts ocean vars (saln,temp,velocity)
! --- at given depth. Data gets dumped in a ascii file.
! --- ---------------------------------------------------
! --- Author: Knut Liseter

program gpdens
use mod_read_gp
use mod_spline_calc
use mod_data_new, only : data_new
use m_density
use m_densitylevels
use m_time_filter
!use m_jot
use m_exceedence
use m_accplot
use m_mkestat
use m_rotate
use m_correlation
!use m_sort
use mod_netcdf_ops
implicit none

#if defined (IARGC) 
integer*4, external :: iargc
#endif

character(len=100) :: infile,outfile,cpost,outfile2, tmparg, stname
character(len= 3) :: rungen
integer :: i,j,ios,rlen,n,irec, ivrec, nrm, ncid
real, allocatable, dimension(:) :: ut,vt,temp,saln,intf
logical :: ex
type(data_new), allocatable :: gp(:), gp2(:), gpfilt(:), gptide(:)
type(data_new), allocatable :: cm(:), cm2(:), cmfilt(:), cmtide(:)
real :: thetagp, maxA, tmp(1), gplon, gplat, cmlon, cmlat, thetacm
real torad,todeg, dnew
integer, parameter :: nlev=15   ! <100
integer, parameter :: maxrec=10000   ! Change if needed
torad=pi/180.0
todeg=1.0/torad



! Get args
if (iargc()==3) then
   call getarg(1,rungen)
   call getarg(2,infile)
   call getarg(3,tmparg) ; read(tmparg,*) dnew
else
   print *,'Routine extracts data from a gp-file for velocity '
   print *,'and computes various statistics for a specified depth'
   print *,'depth level.'
   print *
   print *,'The statistics are dumped to a netcdf file called '
   print *,'gpdens.nc. This file has vertical dimension as unlimited'
   print *,'dimension, meaning that you can add more depth levels '
   print *,'to this file by running gpdens repeatedly. '
   print *
   print *,'Originally this routine was used for comparing model'
   print *,'with observations, currently it only computes model'
   print *,'Statistics. (there are pointers in the code on what'
   print *,'needs to be changed)'
   print *
   print *,'If the routine is run for several times for a infile,'
   print *,'the result will appended  to an already present netcdf file. '
   print *,'****************************************************'
   print *,'Usage: gpdens rungen infile depthlevel'
   print *,'       rungen    : 3-letter ID of model run'
   print *,'       infile    : name of gp-file'
   print *,'       depthlevel: which depth to extract data from'
   stop '(gpseries)'
end if

! Check if file exists
inquire(exist=ex,file=trim(infile))
if (.not. ex) then
   print *,'file '//trim(infile)//' does not exist'
   stop 
end if

! Get info on what is in one gp file record. 
call read_gpheader(rungen)
call  spline_calc_ini_frominput('spline',(/dnew/),1)

! layer values - te
allocate(ut  (kdm))
allocate(vt  (kdm))
allocate(temp(kdm))
allocate(saln(kdm))
allocate(intf(kdm))
allocate(gp (maxrec))

! Read GP records and put into GP variable
ios=0
irec=1  ! index of file records
ivrec=1 ! Index of alid records
do while (ios==0) 
   ! Read data into temporary variable r4s in mod_read_gp
   ! NB: This reads _everything_ in one time record into r4s
   call read_gprecord(infile,irec,ios)
   if (read_success) then

      ! Now Read from the temporary var r4s
      gplon =readvar('lon')
      gplat =readvar('lat')
      ut  =readvar('ut'  ,kdm)
      vt  =readvar('vt'  ,kdm)
      temp=readvar('temp',kdm)
      saln=readvar('saln',kdm)
      intf=readvar('intf',kdm)

      if (intf(kdm)<dnew) then
         print *,'input depth exceeds GP depth'
         print *,intf
         stop
      end if

      ! Verical spline calc
      call spline_calc_1d(ut  ,intf,tmp(1),ndeep,kdm) ;gp(ivrec)%u=tmp(1)
      call spline_calc_1d(vt  ,intf,tmp(1),ndeep,kdm) ;gp(ivrec)%v=tmp(1)
      call spline_calc_1d(temp,intf,tmp(1),ndeep,kdm) ;gp(ivrec)%t=tmp(1)
      call spline_calc_1d(saln,intf,tmp(1),ndeep,kdm) ;gp(ivrec)%s=tmp(1)

      ! m/s -> cm/s
      gp(ivrec)%u=gp(ivrec)%u*100.
      gp(ivrec)%v=gp(ivrec)%v*100.
      gp(ivrec)%speed=sqrt( gp(ivrec)%u**2 + gp(ivrec)%v**2 )
      gp(ivrec)%dir=atan2(gp(ivrec)%u,gp(ivrec)%v)*todeg
      ivrec=ivrec+1
   end if
   irec=irec+1
end do
nrm=ivrec-1
print *,'valid records',nrm

!-----------------> Modify here for comparing model with obs <---------------
!!!!!!!!!!!!!!!!!!!!!!!
! Read observations (cm vars ) here

!!!!!!!!!!!!!!!!!!!!!!!!
! TODO: make sure time variable in gp/cm is correct
!---------------> end modify here for comparing model with obs <---------------


!!!!!!!!!!!!!!!!!!!!!!1
! I have no clue on some of these values...
stname='test'
cmlon=gplon
cmlat=gplat

! Fill netcdf with som initial attributes
i=index(infile,'.uf')
if (i==0) then
   print *,'Error: infile '//trim(infile)//' is not a .uf file?'
   call exit(1)
end if
ncfile=infile(1:i-1)//'_gpdens.nc'
call ncopencreate(ncfile,ncid)
call ncputatt(ncid,'NF90_GLOBAL','GP_longitude',gplon)
call ncputatt(ncid,'NF90_GLOBAL','GP_latitude' ,gplat)
call ncputatt(ncid,'NF90_GLOBAL','CM_longitude',cmlon)
call ncputatt(ncid,'NF90_GLOBAL','CM_latitude' ,cmlat)
call ncputatt(ncid,'NF90_GLOBAL','CM_station_id' ,stname)
call ncerr(NF90_CLOSE(ncid))

!-----------------> Modify here for comparing model with obs <---------------
! TODO : Compute joint timeseries obs/model - needs day to be properly set
!! merging the two data sets
!   call merge(gp,cm,gp2,cm2,ngp,ncm,maxnr,nrm)
! Temoprary cheat for gp2 and cm2
allocate(gp2(nrm))
allocate(cm2(nrm))
allocate(gptide(nrm),gpfilt(nrm))
allocate(cmtide(nrm),cmfilt(nrm))
gp2=gp(1:nrm)
cm2=gp(1:nrm)
print *,nrm
!---------------> end modify here for comparing model with obs <---------------


! generate exceedence diagnostics
call exceedence(cm2,gp2,nrm,nint(dnew))

! generate accumulation plots.
call accplot(cm2,gp2,nrm,nint(dnew))

! compute density functions 
maxA=0.0
call density(cm2,nrm,'CM',thetacm,nint(dnew),maxA,'U')
call density(gp2,nrm,'GP',thetagp,nint(dnew),maxA,'U')

! compute low-pass filtered timeseries
call time_filter(cm2,gp2,cmfilt,gpfilt,cmtide,gptide,nrm,24.0,'AVE',nint(dnew))

! compute density functions for filtered velocities
maxA=0.0
call density(cmfilt,nrm,'CM',thetacm,nint(dnew),maxA,'F')
call density(gpfilt,nrm,'GP',thetagp,nint(dnew),maxA,'F')

! compute density functions for tidal velocities
maxA=0.0
call density(cmtide,nrm,'CM',thetacm,nint(dnew),maxA,'T')
call density(gptide,nrm,'GP',thetagp,nint(dnew),maxA,'T')

! compute MKE statistics
call mkestat(cm2,gp2,cmfilt,gpfilt,cmtide,gptide,nrm,nint(dnew)) !  test

! complex correlation
call correlation(cm2,gp2,cmfilt,gpfilt,cmtide,gptide,nrm,nint(dnew)) ! test

! print low pass filtered timeseries rotated to principal direction
call rotate(cmfilt,gpfilt,nrm,thetacm,thetagp,nint(dnew)) ! test

print *,'Normal exit: data in '//trim(ncfile)

end program

