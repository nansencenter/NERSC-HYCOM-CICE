! --- ---------------------------------------------------
! --- Simple program for reading generic gp series
! --- and preparing them for tecplot. Produces
! --- a time/vs depth (hovmoller) plot best viewed with
! --- contour  in tecplot. No vertical interpolation!
! ---
! --- Files needed to run this program are:
! ---  1) a gp-unformatted file (created by HYCOM)
! ---  2) a gp header file      (created by HYCOM)
! ---
! --- Header file contains info on what data can be found 
! --- in one gp-record. This makes this routine 
! --- independent of model size erc...
! ---
! --- ---------------------------------------------------
! --- Author: Knut Liseter

program gphovmoller
use mod_read_gp
use mod_year_info
!use m_masslayers
use netcdf
use m_ncerr
implicit none
#if defined (IARGC) 
integer*4, external :: iargc
#endif
real, parameter :: undef=-999.
character(len=80) :: tmparg
character(len=200) :: infile,outfile,outfilenc,gpdir, infile2
character(len= 3) :: rungen
character(len= 4) :: cdd
integer :: i,j,ios,rlen,n,n2, k,ind
logical :: ex
real   , allocatable, dimension(:) :: intf, temp, saln, ut, vt, spd
integer :: nip
integer :: i_temp, i_saln,i_intf,i_ut,i_vt
logical :: allequal
integer :: diy, diy_now
integer :: irec, irec2
integer :: ncid, dimid, idtime, idyear, idmonth, iddom, idhour,idtemp, &
           idsaln, idut, idvt, idspd, kdmid, idintf


! Get args
if (iargc()==3) then
   call getarg(1,rungen)
   call getarg(2,infile)
   call getarg(3,gpdir)
else
   print *,'Routine exctract hovmoller data for the GP station '
   print *,'and dumps it in a netcdf file. Variables extracted '
   print *,'are temperature, salinity, current, lower layer '
   print *,'interface and current speed. No vertical interpolation '
   print *,'is done here.'
   print *
   print *,'Usage: gphovmoller rungen infile directory'
   print *,'       rungen: 3-letter ID of model run'
   print *,'       infile: name of gp-file'
   print *,'       directory: location of gp-file'
   stop '(gpseries)'
end if

! Get info on what is in one gp file record. 
infile2=trim(gpdir)//'/'//trim(infile)
call read_gpheader(rungen,trim(gpdir))

! Allocata data holders
allocate(ut  (kdm))
allocate(vt  (kdm))
allocate(temp(kdm))
allocate(saln(kdm))
allocate(intf(kdm))
allocate(spd (kdm))


! Open outfile (based on infile name)
outfile=infile
outfilenc=infile
ind=index(outfile,'.uf')
outfile(ind:ind+4) = '.tec'
outfilenc(ind:ind+3) = '.nc'
outfilenc = 'hovmoller'//trim(outfilenc)

! Define netcdf dims and vars
call ncerr(NF90_create(trim(outfilenc),NF90_CLOBBER,ncid))
call ncerr(NF90_def_dim(ncid,'time',NF90_UNLIMITED,dimid))
call ncerr(NF90_def_dim(ncid,'kdm',kdm,kdmid))
call ncerr(NF90_def_var(ncid,'time',NF90_FLOAT,(/dimid/),idtime))
call ncerr(NF90_def_var(ncid,'year',NF90_INT,(/dimid/),idyear))
call ncerr(NF90_def_var(ncid,'month',NF90_INT,(/dimid/),idmonth))
call ncerr(NF90_def_var(ncid,'day_of_month',NF90_INT,(/dimid/),iddom))
call ncerr(NF90_def_var(ncid,'hour',NF90_INT,(/dimid/),idhour))
call ncerr(NF90_def_var(ncid,'intf',NF90_FLOAt,(/kdmid,dimid/),idintf))
call ncerr(NF90_def_var(ncid,'temp',NF90_FLOAt,(/kdmid,dimid/),idtemp))
call ncerr(NF90_def_var(ncid,'saln',NF90_FLOAt,(/kdmid,dimid/),idsaln))
call ncerr(NF90_def_var(ncid,'ut',NF90_FLOAt,(/kdmid,dimid/),idut))
call ncerr(NF90_def_var(ncid,'vt',NF90_FLOAt,(/kdmid,dimid/),idvt))
call ncerr(NF90_def_var(ncid,'speed',NF90_FLOAt,(/kdmid,dimid/),idspd))
call ncerr(NF90_enddef(ncid))

irec=1
irec2=1
ios=0
do while(ios==0)
   ! Read into temporary array r4s in mod_read_gp
   call read_gprecord(trim(infile2),irec,ios)
   if (ios==0.and.itime>0) then

      ! Get 3d vars from temporary array r4s in mod_read_gp
      intf(:) = readvar('intf',kdm)
      temp(:) = readvar('temp',kdm)
      saln(:) = readvar('saln',kdm)
      ut  (:) = readvar('ut'  ,kdm)
      vt  (:) = readvar('vt'  ,kdm)
      spd=sqrt(ut**2 +vt**2)

      call ncerr(NF90_put_var(ncid,idtime ,fyear ,start=(/1,irec2/)))
      call ncerr(NF90_put_var(ncid,idyear ,iyear ,start=(/1,irec2/)))
      call ncerr(NF90_put_var(ncid,idmonth,imonth,start=(/1,irec2/)))
      call ncerr(NF90_put_var(ncid,iddom  ,iday+1,start=(/1,irec2/)))
      call ncerr(NF90_put_var(ncid,idhour ,ihour ,start=(/1,irec2/)))
      call ncerr(NF90_put_var(ncid,idintf ,intf  ,start=(/1,irec2/)))
      call ncerr(NF90_put_var(ncid,idtemp ,temp  ,start=(/1,irec2/)))
      call ncerr(NF90_put_var(ncid,idsaln ,saln  ,start=(/1,irec2/)))
      call ncerr(NF90_put_var(ncid,idut   ,ut    ,start=(/1,irec2/)))
      call ncerr(NF90_put_var(ncid,idvt   ,vt    ,start=(/1,irec2/)))
      call ncerr(NF90_put_var(ncid,idspd  ,spd   ,start=(/1,irec2/)))

      irec2=irec2+1
   end if
   irec=irec+1
end do
call ncerr(NF90_close(ncid))
print *,'Diag dumped to netcdf file '//trim(outfilenc)
end program

