      program p_mldcompare
      !
      ! This program must be linked against netcdf
      use mod_xc, only: idm,jdm, xcspmd
      use mod_za, only: zaiost
      use m_dbm_mld_read
      use mod_hycomfile_io
      use mod_sphere_tools
      use mod_grid 
      use mod_confmap
      use mod_year_info
      use netcdf
      use mod_regions
      implicit none

      real, parameter :: undef2=-1e14
#if defined (IARGC)
   integer*4, external :: iargc
#endif

      character(len=80) :: hycfile, fbase
      integer :: nrmem

      integer       :: nblon, nblat
      real, pointer :: obsmld(:,:),mldlon(:),mldlat(:)


      character(len=3) :: creg
      integer, allocatable :: regflag(:,:)
      real, allocatable :: modmld(:,:), temp(:,:), dp(:,:), dpsum(:,:), nmodmld(:,:) , temp1(:,:)
      character(len=11) tag7
      logical :: ex,regionflag,havedata
      integer :: i,j,ipiv,jpiv,jday,ireg,diy
      real :: lat_n,lon_n


      integer :: stat,k
      integer :: ncid, dimid_idm, dimid_jdm, obs_varid, mod_varid, maskvarid
      logical, allocatable :: mask(:,:)
      real :: rdummy
      type(hycomfile) :: hfile


      ! Regional stats
      integer :: ios, iosa, iosb, iosc, ierr, kdm
      character(len=80) :: regname, filetype
      integer :: numregion
      integer :: npoints, npoints2
      real :: diffsq, meanobs, meanmod, RMS, bias, yfrac


      ! Process input args
      if (iargc()/=1) then
          print *,'Usage:'
          print *,'  slacmp hycom-file'
          stop
       else
          call getarg(1,hycfile)
       end if

      ! Deduce hycom filetype
      filetype= getfiletype(hycfile)

      ! Initialize the hycomfile module - and retrieve vertical levels
      call initHF(hfile,trim(hycfile),trim(filetype))
      kdm=vdim(hfile)

      ! get julian day rel this year-1-1
      jday=hfile%iday
      diy =datetojulian(hfile%iyear,12,31,hfile%iyear,1,1)

      ! "Float year"
      yfrac=hfile%iyear + float(jday)/float(diy)

      print *,'Julian day is ',jday



      ! Read duacs data from opendap
      call dBM_MLD_read(obsmld,mldlon,mldlat,nblon,nblat,undef2,jday)
      print *

      ! Get model grid
      call xcspmd()
      call zaiost()
      call get_grid() ; where(depths>0.5*huge) depths=0.
      call initconfmap(idm,jdm)

      ! Read model temp field
      allocate(nmodmld(nblon,nblat))
      allocate(modmld(idm,jdm))
      allocate(temp  (idm,jdm))
      allocate(temp1 (idm,jdm))
      allocate(dp    (idm,jdm))
      allocate(dpsum (idm,jdm))

      ! A bit complicated because we must cycle through layers to get MLD
      k=1
      modmld=0.
      temp=0.
      dpsum=0.
      do while (k<=kdm)

         call HFReadDPField(hfile,dp,idm,jdm,k,1)
         call HFReadField(hfile,temp,idm,jdm,'temp    ',k,1)
         dp=dp/onem

         if (k==1) temp1=temp
         dpsum=dpsum+dp

         where(depths>.1 .and. abs(temp-temp1)>0.2 .and. modmld==0.)
            modmld=dpsum
         end where
         k=k+1
      end do

      where(depths>.1 .and. modmld==0.)
         modmld=dpsum
      end where


         
        




      ! Find model points corr to mld points
      nmodmld=undef2
      do i=1,nblon
      do j=1,nblat
         call oldtonew(mldlat(j),mldlon(i),lat_n,lon_n)
         call pivotp(lon_n,lat_n,ipiv,jpiv)
         if ( ipiv>1 .and. ipiv<idm .and.  jpiv>1 .and. jpiv<jdm ) then
            if (depths(ipiv,jpiv)>1. ) then
               nmodmld(i,j)=modmld(ipiv,jpiv)   
            end if
         end if
      end do
      end do

      ! Correction factor for pressure coords
      if (maxval(nmodmld,mask=nmodmld>undef2) > 10.*9806) then
         where(nmodmld>undef2) nmodmld=nmodmld/9806.
      end if


      ! Quick dump of data points to netcdf file
      stat = NF90_CREATE('mldcmp.nc', NF90_CLOBBER , ncid)
      if (stat /= nf90_noerr) then
         print *,nf90_strerror(stat)
         call exit(3)
      end if
      stat = NF90_DEF_DIM (NCID, 'nblon', nblon, dimid_idm)
      stat = NF90_DEF_DIM (NCID, 'nblat', nblat, dimid_jdm)
      stat = NF90_DEF_VAR(NCID, 'mldobs', NF90_Float, (/dimid_idm,dimid_jdm/), obs_varid)
      stat = NF90_DEF_VAR(NCID, 'mldmod', NF90_Float, (/dimid_idm,dimid_jdm/), mod_varid)
      stat = NF90_PUT_ATT(NCID, obs_varid, '_FillValue', real(undef2,kind=4) )
      stat = NF90_PUT_ATT(NCID, mod_varid, '_FillValue', real(undef2,kind=4) )
      STAT = NF90_ENDDEF(NCID)
      stat = NF90_PUT_VAR(NCID,mod_VARID,nmodmld,(/1,1/),(/nblon,nblat/))
      if (stat /= nf90_noerr) then
         print *,nf90_strerror(stat)
         call exit(3)
      end if
      stat = NF90_PUT_VAR(NCID,obs_varid,obsmld,(/1,1/),(/nblon,nblat/))
      if (stat /= nf90_noerr) then
         print *,nf90_strerror(stat)
         call exit(3)
      end if


      ! Initialize regions for RMS calc
      call read_regions()
      numregion=getnregions()
      print *,numregion
      allocate(mask   (nblon,nblat))
      allocate(regflag(nblon,nblat))
      do ireg=1,numregion

         print *,'iregion: ',ireg

         ! Get mask for this region
         call getmask(ireg,mask,mldlon,mldlat,nblon,nblat)

         where (mask .and. obsmld/=undef2 .and. nmodmld/=undef2)
            mask=.true.
         elsewhere
            mask=.false.
         end where
         print *,count(mask)

         meanmod=0.
         meanobs=0.
         diffsq=0.
         npoints=0
         do j=1,nblat
         do i=1,nblon
         if (mask(i,j)) then
            npoints=npoints+1
            meanobs=meanobs+obsmld(i,j)
            meanmod=meanmod+nmodmld(i,j)
            diffsq=diffsq+(obsmld(i,j)-nmodmld(i,j))**2
            regflag(i,j)=regflag(i,j)+2**ireg ! Can be binary decoded
         end if
         end do
         end do

         bias=(meanmod-meanobs)/npoints
         RMS = sqrt( diffsq/npoints - bias**2)
         meanmod=meanmod/npoints
         meanobs=meanobs/npoints
         call dumprms_ascii(ireg,yfrac,meanmod,meanobs,bias,rms,npoints,'mldRMS_')
         print *,'meanmod:',meanmod
         print *,'meanobs:',meanobs
         print *,'bias   :',bias   
         print *,'RMS    :',RMS    
      end do


      if (numregion>0) then

         ! put mask in netcdf file
         STAT = NF90_REDEF(NCID)
         stat = NF90_DEF_VAR(NCID, 'regionmask',NF90_INT,(/dimid_idm,dimid_jdm/), maskvarid)
         if (stat /= nf90_noerr) then
            print *,nf90_strerror(stat)
            call exit(3)
         end if
         STAT = NF90_ENDDEF(NCID)

         stat = NF90_PUT_VAR(NCID,maskvarid,regflag,(/1,1/),(/nblon,nblat/))
         if (stat /= nf90_noerr) then
            print *,nf90_strerror(stat)
            call exit(3)
         end if
      end if


      stat = NF90_CLOSE(ncid)


      print *
      print *,'Data dumped to:'
      do ireg=1,numregion
         print *,trim('mldRMS_')//trim(getname(ireg))//'.asc'
      end do
      print *,'mldcmp.nc'



      end program

