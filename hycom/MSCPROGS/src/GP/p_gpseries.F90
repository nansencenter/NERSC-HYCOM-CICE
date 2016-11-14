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

program gpseries
use mod_read_gp
use mod_year_info
!use m_calcconst
!use m_calcu
use mod_spline_calc
use m_masslayers
use netcdf
use m_ncerr
implicit none
#if defined (IARGC)
integer*4, external :: iargc
#endif
real, parameter :: undef999=-999.
character(len=80) :: tmparg
character(len=80) :: infile,outfile, ncfile,gpdir, infile2
character(len= 3) :: rungen
character(len= 4) :: cdd
integer :: i,j,ios,rlen,n
real, allocatable, dimension(:) :: a,b,c
integer :: k,ind
integer :: ic
logical :: ex
real    :: dnew,unew,vnew,snew,tnew, spd, tmp(1)
integer, allocatable, dimension(:) :: ipp
real   , allocatable, dimension(:) :: intf, temp, saln, ut, vt
integer :: nip
logical :: allequal
integer :: diy, diy_now
integer :: irec, irec2
integer :: ncid, dimid, idtime, idyear, idmonth, iddom, idhour,idtemp, &
           idsaln, idut, idvt, idspd


! Get args
if (iargc()==4) then
   call getarg(1,rungen)
   call getarg(2,infile)
   call getarg(3,tmparg) ; read(tmparg,*) dnew
   call getarg(4,gpdir) ;
else
   print *,'Extracts pre-set variables from a GP file and dumps them in netcdf format.'
   print *,'It extracts temperature, salinity, currents and current speed for the '
   print *,'GP station specified'
   print *
   print *,'Usage: gpseries rungen infile depthlevel [directory]'
   print *,'       rungen    : 3-letter ID of model run'
   print *,'       infile    : name of gp-file'
   print *,'       depthlevel: which depth to extract data from'
   print *,'       directory : where gp data is located'
   stop '(gpseries)'
end if

! Check if file exists
infile2=trim(gpdir)//'/'//trim(infile)
inquire(exist=ex,file=infile2)
if (.not. ex) then
   print *,'file '//trim(infile2)//' does not exist'
   stop 
end if

! Get info on what is in one gp file record. 
call read_gpheader(rungen,trim(gpdir))
print *,'kdm=',kdm
allocate(ut  (kdm))
allocate(vt  (kdm))
allocate(temp(kdm))
allocate(saln(kdm))
allocate(intf(kdm))
!END!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Spline constants
allocate(a(kdm))
allocate(b(kdm))
allocate(c(0:kdm+1))
allocate(ipp(kdm))

! Open outfile (based on infile name)
outfile=infile
ind=index(outfile,'.uf')
write(cdd,'(i4.4)') int(dnew)
outfile(ind:ind+9) = '_d'//cdd//'.asc'
print *
open(11,file=trim(outfile),status='unknown',form='formatted')

! New - put results into netcdf file
ncfile=infile
ncfile(ind:ind+8) = '_d'//cdd//'.nc'
call ncerr(NF90_create(trim(ncfile),NF90_CLOBBER,ncid))
call ncerr(NF90_def_dim(ncid,'time',NF90_UNLIMITED,dimid))
call ncerr(NF90_def_var(ncid,'time',NF90_FLOAT,(/dimid/),idtime))
call ncerr(NF90_def_var(ncid,'year',NF90_INT,(/dimid/),idyear))
call ncerr(NF90_def_var(ncid,'month',NF90_INT,(/dimid/),idmonth))
call ncerr(NF90_def_var(ncid,'day_of_month',NF90_INT,(/dimid/),iddom))
call ncerr(NF90_def_var(ncid,'hour',NF90_INT,(/dimid/),idhour))
call ncerr(NF90_def_var(ncid,'temp',NF90_FLOAt,(/dimid/),idtemp))
call ncerr(NF90_def_var(ncid,'saln',NF90_FLOAt,(/dimid/),idsaln))
call ncerr(NF90_def_var(ncid,'ut',NF90_FLOAt,(/dimid/),idut))
call ncerr(NF90_def_var(ncid,'vt',NF90_FLOAt,(/dimid/),idvt))
call ncerr(NF90_def_var(ncid,'speed',NF90_FLOAt,(/dimid/),idspd))
call ncerr(NF90_enddef(ncid))

call  spline_calc_ini_frominput('spline',(/dnew/),1)


n=1
ios=0
irec=1
irec2=1 ! Number of non-empty records
do while(ios==0)

   !  Read gprecord into temporary array r4s in mod_read_gp
   call read_gprecord(infile2,irec,ios)
   if (ios==0.and.itime>0) then

      ! Get 3d vars from temporary array r4s in mod_read_gp
      ut  =readvar1d('ut'  ,kdm)
      vt  =readvar1d('vt'  ,kdm)
      temp=readvar1d('temp',kdm)
      saln=readvar1d('saln',kdm)
      intf=readvar1d('intf',kdm)

      ! ipp points to mass-filled layers
      call masslayers(intf,kdm,ipp,nip)

      do k=1,nip
         intf(k)=intf(ipp(k))
         temp(k)=temp(ipp(k))
         saln(k)=saln(ipp(k))
         ut  (k)=ut  (ipp(k))
         vt  (k)=vt  (ipp(k))
      end do

!      call calcconst(nip,kdm,ut,intf,a,b,c,'A','a')
!      unew=calcu(nip,kdm,a,b,c,intf,dnew,ic,undef999)  
!
!      call calcconst(nip,kdm,vt,intf,a,b,c,'A','a')
!      vnew=calcu(nip,kdm,a,b,c,intf,dnew,ic,undef999)  
!
!      call calcconst(nip,kdm,temp,intf,a,b,c,'A','a')
!      tnew=calcu(nip,kdm,a,b,c,intf,dnew,ic,undef999)  
!
!      call calcconst(nip,kdm,saln,intf,a,b,c,'A','a')
!      snew=calcu(nip,kdm,a,b,c,intf,dnew,ic,undef999)

      call spline_calc_1d(ut(1:nip)  ,intf(1:nip),tmp(1),ndeep,nip); unew=tmp(1)
      call spline_calc_1d(vt(1:nip)  ,intf(1:nip),tmp(1),ndeep,nip); vnew=tmp(1)
      call spline_calc_1d(temp(1:nip),intf(1:nip),tmp(1),ndeep,nip); tnew=tmp(1)
      call spline_calc_1d(saln(1:nip),intf(1:nip),tmp(1),ndeep,nip); snew=tmp(1)


      spd=sqrt(unew**2+vnew**2)

      write(11,'(f14.6,4i5," ",6f10.4)')  &
         fyear,iyear,imonth,iday,ihour,dnew,snew,tnew,unew,vnew

      call ncerr(NF90_put_var(ncid,idtime ,fyear ,start=(/irec2/)))
      call ncerr(NF90_put_var(ncid,idyear ,iyear ,start=(/irec2/)))
      call ncerr(NF90_put_var(ncid,idmonth,imonth,start=(/irec2/)))
      call ncerr(NF90_put_var(ncid,iddom  ,iday+1,start=(/irec2/)))
      call ncerr(NF90_put_var(ncid,idhour ,ihour ,start=(/irec2/)))
      call ncerr(NF90_put_var(ncid,idtemp ,tnew  ,start=(/irec2/)))
      call ncerr(NF90_put_var(ncid,idsaln ,snew  ,start=(/irec2/)))
      call ncerr(NF90_put_var(ncid,idut   ,unew  ,start=(/irec2/)))
      call ncerr(NF90_put_var(ncid,idvt   ,vnew  ,start=(/irec2/)))
      call ncerr(NF90_put_var(ncid,idspd  ,spd   ,start=(/irec2/)))

      irec2=irec2+1
   end if
   irec=irec+1
end do
close(10)
close(11)
call ncerr(NF90_close(ncid))

! Copy to a generic file as well, useful for scripts.
call system("cp "//trim(outfile)//" gp.asc")
print '(a)','Results dumped to '//trim(outfile)//", (copied to gp.asc as well)"
print '(a)','Results dumped to netcdf file '//trim(ncfile)

end program

