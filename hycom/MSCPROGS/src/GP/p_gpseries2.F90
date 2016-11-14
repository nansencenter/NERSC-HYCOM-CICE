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

program gpseries2
use mod_read_gp
use mod_year_info
use mod_spline_calc
!use m_calcconst
!use m_calcu
use m_masslayers
use m_ncerr
use netcdf
implicit none
#if defined (IARGC) 
integer*4, external :: iargc
#endif
real, parameter :: undef999=-999.
character(len=80) :: tmparg
character(len=100) :: infile,outfile,cpost,outfile2,ncfile,infile2,gpdir
character(len=20) :: cvar
character(len= 3) :: rungen
character(len= 4) :: cdd
integer :: i,j,ios,rlen,n
real, allocatable, dimension(:) :: a,b,c
integer :: k,ind
integer :: ic
logical :: ex
real    :: dnew,unew,vnew,snew,tnew, tmp(1)
integer, allocatable, dimension(:) :: ipp
real   , allocatable, dimension(:) :: var,intf
integer :: nip
integer :: i_temp, i_saln,i_intf,i_ut,i_vt,i_cvar
real    :: newvar
logical :: allequal,ocnvar
integer :: irec, irec2
integer :: ncid, dimid, idtime, idyear, idmonth, iddom, idhour,idtemp, &
           idsaln, idut, idvt, idspd


! Get args
if (iargc()==5) then
   call getarg(1,rungen)
   call getarg(2,infile)
   call getarg(3,cvar)
   call getarg(4,tmparg) ; read(tmparg,*) dnew
   call getarg(5,gpdir)
else
   print *,'Routine extracts data from a gp-file for the given '
   print *,'variable name and given gp-file. Can extract both  '
   print *,'2D vars and 3D vars (depthlevel only applicable    '
   print *,'to 3D vars....'
   print *,'****************************************************'
   print *,'Usage: gpseries2 rungen infile variable depthlevel '
   print *,'       rungen    : 3-letter ID of model run'
   print *,'       infile    : name of gp-file'
   print *,'       variable  : name of variable to extract'
   print *,'       depthlevel: which depth to extract data from (ignored for 2D vars)'
   print *,'       directory : Directory containing gp data  '
   print *
   stop '(gpseries)'
end if

! Check if file exists
infile2=trim(gpdir)//'/'//trim(infile)
inquire(exist=ex,file=trim(infile2))
if (.not. ex) then
   print *,'file '//trim(infile2)//' does not exist'
   stop 
end if

! Get info on what is in one gp file record. 
call read_gpheader(rungen,trim(gpdir))

! Get  variable
!BEGIN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
i=getvarind('intf')
if (i==-1) then
   print *,'Could not get variable intf'
   stop
else
   allocate(intf(vardim(i)))
   i_intf=i
end if
kdm=vardim(i_intf)

i=getvarind(trim(cvar))
if (i==-1) then
   print *,'Could not get variable '//trim(cvar)
   stop
else
   allocate(var(vardim(i)))
   i_cvar=i
end if
!END!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ocnvar=vardim(i_intf)==vardim(i_cvar)
print *
if (ocnvar) then
   print *,'This is a 3D variable-extracting at depth ',dnew
else if (vardim(i_cvar)==1) then
   dnew=0.
   if (iargc()==4) then
      write(6,*) 'This is a 2D variable  - depth set to zero', &
              '---->supplied depth level was ignored'
   else
      write(6,*) 'This is a 2D variable  - depth set to zero'
   end if
else 
   print *,'This is a 3D variable, but dim is different from ocean ...'
   stop '(gpseries2)'
end if
print *


! Spline constants
allocate(ipp(kdm))
allocate(a(kdm))
allocate(b(kdm))
allocate(c(0:kdm+1))


! Open outfile (based on infile name)
write(cdd,'(i4.4)') int(dnew)
cpost='_'//trim(cvar)//'_d'//cdd//'.asc'

outfile=infile
ind=index(outfile,'.uf')
outfile(ind:ind+len_trim(cpost))=trim(cpost)
open(11,file=trim(outfile),status='unknown',form='formatted')

! New - put results into netcdf file
ncfile=infile
cpost='_'//trim(cvar)//'_d'//cdd//'.nc'
ncfile(ind:ind+len_trim(cpost))=trim(cpost)
call ncerr(NF90_create(trim(ncfile),NF90_CLOBBER,ncid))
call ncerr(NF90_def_dim(ncid,'time',NF90_UNLIMITED,dimid))
call ncerr(NF90_def_var(ncid,'time',NF90_FLOAT,(/dimid/),idtime))
call ncerr(NF90_def_var(ncid,'year',NF90_INT,(/dimid/),idyear))
call ncerr(NF90_def_var(ncid,'month',NF90_INT,(/dimid/),idmonth))
call ncerr(NF90_def_var(ncid,'day_of_month',NF90_INT,(/dimid/),iddom))
call ncerr(NF90_def_var(ncid,'hour',NF90_INT,(/dimid/),idhour))
call ncerr(NF90_def_var(ncid,trim(cvar),NF90_FLOAT,(/dimid/),idtemp))
call ncerr(NF90_enddef(ncid))

call  spline_calc_ini_frominput('spline',(/dnew/),1)

ios=0
irec=1
irec2=1
do while(ios==0)
   call read_gprecord(infile2,irec,ios)
   !print *,' itime tst:',itime
   if (ios==0.and.itime>0) then


      call ncerr(NF90_put_var(ncid,idtime ,fyear ,start=(/irec2/)))
      call ncerr(NF90_put_var(ncid,idyear ,iyear ,start=(/irec2/)))
      call ncerr(NF90_put_var(ncid,idmonth,imonth,start=(/irec2/)))
      call ncerr(NF90_put_var(ncid,iddom  ,iday+1,start=(/irec2/)))
      call ncerr(NF90_put_var(ncid,idhour ,ihour ,start=(/irec2/)))
      call ncerr(NF90_put_var(ncid,idtemp ,tnew  ,start=(/irec2/)))

      ! Get 3d vars
      if (ocnvar) then
         var  =readvar(cvar  ,kdm)
         intf =readvar('intf',kdm)

         ! ipp points to mass-filled layers
         call masslayers(intf,kdm,ipp,nip)

         do k=1,nip
            intf(k)=intf(ipp(k))
            var (k)=var (ipp(k))
         end do

         !call calcconst(nip,kdm,var,intf,a,b,c,'A','a')
         !newvar=calcu(nip,kdm,a,b,c,intf,dnew,ic,undef999)  
         call spline_calc_1d(var(1:nip)  ,intf(1:nip),tmp(1),ndeep,nip); newvar=tmp(1)

         write(11,'(f14.6,4i5," ",2f10.4)')  &
            fyear,iyear,imonth,iday,ihour,dnew,newvar
         call ncerr(NF90_put_var(ncid,idtemp ,newvar  ,start=(/irec2/)))
      else if (vardim(i_cvar)==1) then
         var = readvar(cvar)
         write(11,'(f14.6,4i5," ",2f10.4)')  &
            fyear,iyear,imonth,iday,ihour,dnew,var(:)
         call ncerr(NF90_put_var(ncid,idtemp ,var  ,start=(/irec2/)))
      end if
      irec2=irec2+1
   end if
   irec=irec+1
end do
close(10)
close(11)
call ncerr(nf90_close(ncid))

! Copy to a generic file as well, useful for scripts.
outfile2='gp_'//trim(cvar)//'_d'//trim(cdd)//'.asc'
call system("cp "//trim(outfile)//" "//trim(outfile2))

print '(a)','Results dumped to '//trim(outfile)//", also copied to ", &
            trim(outfile2)
print '(a)','Results dumped to netcdf file '//trim(ncfile)


end program

