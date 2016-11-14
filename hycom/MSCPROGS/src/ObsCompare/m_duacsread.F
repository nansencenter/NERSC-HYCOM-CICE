      module m_duacsread
      contains
      subroutine duacsread(slaout,lonout,latout,nlonin,nlatin,undef,
     &   jday,fnamein)
      ! Input arrays are hardcoded for now - pointer is a viable option
      !
      ! Routine to compare modeled sla with observations 
      ! Program takes three arguments  - first is model
      ! field (only DODS link for now), second is SLA
      ! field - also DODS. Third is some form of time 
      ! unit (to be decided)
      !
      ! This program must be linked against nc-dap
      use mod_year_info
      implicit none
      !include 'nc-dap.inc'
      include 'netcdf.inc'

      integer, intent(in) :: nlonin, nlatin, jday
      real   , intent(out) :: slaout(nlatin,nlonin)
      real   , intent(out) :: lonout(nlonin)
      real   , intent(out) :: latout(nlatin)
      real   , intent(in)  :: undef
      character(len=*), optional   , intent(in)  :: fnamein

C     character(len=*),parameter ::
C    &   fnameod   ="http://opendap.aviso.oceanobs.com"//
C    &   "/thredds/dodsC/dt_ref_global_merged_msla_h"
C     character(len=*),parameter ::
C    &   fnameod   ="http://opendap.aviso.oceanobs.com"//
C    &   "/thredds/dodsC/duacs_global_nrt_msla_merged_h"
C --- OpenDAP path
      character(len=*),parameter :: 
     &   fnameod= "http://opendap.aviso.oceanobs.com/thredds/"//
     &          "dodsC/dt_ref_global_merged_msla_h"
      character(len=*),parameter :: vname="Grid_0001"

      integer :: ncid, stat, ndims, nvars, ngatts,varid,dvarid
      integer :: unlimdimid
      character(len=NF_MAX_NAME) :: varname
      integer, dimension(NF_MAX_VAR_DIMS) :: vardimids, vardimlen,
     &   varcoordtype
      character(len=NF_MAX_NAME), dimension(nf_max_var_dims) :: 
     &   vardimname
      integer :: xtype,varndims,nvaratts
      integer :: i,dimid(1),j,k,mintdiff,tmp
      double precision,  allocatable :: lats(:),lons(:)
      real*4,  allocatable :: sla(:,:)
      integer,  allocatable :: times(:)
      integer :: nlons,nlats, ntimes
      integer, allocatable :: indexes(:)
      integer :: timeind,start(3),count(3)
      integer :: dimlen, year, month, day
      real :: fillval
      character(len=200) :: fname

      ! Use duacs opendap in this case
      if (.not.present(fnamein)) then
         fname=trim(fnameod)
      else
         fname=trim(fnamein)
      end if

      stat = NF_OPEN(trim(fname), NF_NOWRITE , ncid)
      if (stat /= nf_noerr) then
         print *,'Error Reading '//trim(fname)
         print *,nf_strerror(stat)
         print *,'This could happen if you try OPENDAP approach'
         print *,'without actual opendap support in netcdf libraries'
         call exit(3)
      end if
      print *,'Reading SLA from ',trim(fname)


      ! Get coordinate variable lat
      stat = NF_INQ_VARID (NCID, "NbLatitudes", varid)
      stat = NF_INQ_dimID (NCID, "NbLatitudes", dimid)
      stat = NF_INQ_dimlen(NCID, dimid, dimlen)
      stat = NF_INQ_dimlen(NCID, dimid, dimlen)
      if (stat /= nf_noerr) then
         print *,nf_strerror(stat)
         call exit(3)
      end if
      if (dimlen == nlatin) then
         allocate(lats    (nlatin))
         stat =  NF_get_VAR_double(NCID, VARID, lats)
         latout=lats
      else
         print *,'Latitude dimension mismatch'
         call exit(2)
      end if

      ! Get coordinate variable lon
      stat = NF_INQ_VARID (NCID, "NbLongitudes", varid)
      stat = NF_INQ_dimID (NCID, "NbLongitudes", dimid)
      stat = NF_INQ_dimlen(NCID, dimid, dimlen)
      if (stat /= nf_noerr) then
         print *,nf_strerror(stat)
         call exit(3)
      end if
      if (dimlen == nlonin) then
         allocate(lons    (nlonin))
         stat =  NF_get_VAR_double(NCID, VARID, lons)
         lonout=lons
      else
         print *,'Latitude dimension mismatch'
         call exit(2)
      end if


      ! NB - Netcdf Files have no time dimension 
      if (.not.present(fnamein)) then
         ! Get coordinate variable time
         stat = NF_INQ_VARID (NCID, "time", varid)
         stat = NF_INQ_dimID (NCID, "time", dimid)
         stat = NF_INQ_dimlen(NCID, dimid, dimlen)
         if (stat /= nf_noerr) then
            print *,nf_strerror(stat)
            call exit(3)
         end if
         ntimes=dimlen
         allocate(times    (dimlen))
         stat =  NF_get_VAR_int(NCID, VARID, times)


         ! Find time index corresponding to julian date jday
         mintdiff=1000000
         do k=1,ntimes
            tmp=abs(jday-times(k))
            if (tmp < mintdiff ) then
               mintdiff=tmp
               timeind=k
            end if
         end do
         if (mintdiff>10) then
            print *,'SLA DUACS: found no times closer than  10 days'
            call juliantodate(times(timeind),year,month,day,
     &                        1950,1,1)
            print *,'Closest match is ',mintdiff, ' at ',year, month,day
            call exit(1)
         else
            call juliantodate(times(timeind),year,month,day,
     &                        1950,1,1)
            print *,'Closest match is ',mintdiff, ' at ',year, month,day
         end if
      end if

         
      ! Testing
      allocate(sla(nlatin,nlonin))
      if (present(fnamein)) then
         count=(/nlatin,nlonin,0/)
         start=(/1,1,0/)
         stat = NF_INQ_VARID (NCID, "Grid_0001", varid)
         stat = NF_GET_VARA_REAL(NCID,VARID,START(1:2),COUNT(1:2),sla)
      else
         count=(/nlatin,nlonin,1/)
         start=(/1,1,timeind/)
         stat = NF_INQ_VARID (NCID, "Grid_0001", varid)
         stat = NF_GET_VARA_REAL(NCID, VARID, START, COUNT, sla)
      end if

      if (stat /= nf_noerr) then
         print *,nf_strerror(stat)
         call exit(3)
      end if

      stat =  NF_GET_ATT_REAL   (NCID, VARID, '_FillValue',fillval)
      if (stat /= nf_noerr) then
         print *,nf_strerror(stat)
         call exit(4)
      end if


      slaout=sla
      where (sla==fillval) slaout=undef
      if (present(fnamein)) then
         deallocate(lons,lats,sla)
      else
         deallocate(times,lons,lats,sla)
      end if

c --- Hmmm .. Ice values seem to be out of undef range
      where (slaout > 1e18) slaout=undef


      end subroutine
      end module

