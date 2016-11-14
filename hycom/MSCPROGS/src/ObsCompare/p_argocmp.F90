program p_argocmp
   use mod_year_info
   use mod_levitus
   use mod_confmap
   use mod_hycomfile_io
   use netcdf
   use mod_grid
   use mod_xc
   use mod_za 
   use mod_spline_calc
   use m_read_coriolis
   use m_dump_coriolis
   implicit none

#if defined (IARGC)
   integer*4, external :: iargc
#endif

   ! Version info of this program
   character(len=*), parameter :: cver = 'V0.2'

   ! Input file with days
   character(len=80) :: hyc_file
   integer            :: year,day,hour,find,ios,ifile, &
                         counter, refyear
   logical            :: ex

   real, allocatable, dimension(:,:  ) :: dp, pres, presold, saln, temp
   logical, allocatable, dimension(:) :: valid_profile

   character(len=8), dimension(:), pointer   :: floatid
   character(len=8)   :: cdate
   character(len=80)  :: modcorfile,corfile    ! Name of new netcdf file
   character(len=80)  :: cprof
   character(len=40)  :: tmpchar,vartitle,fname
   character(len=  4) :: char_fac
   character(len=3) :: cgg,cii
   character(len=12) :: clon, clat, cjuld
   character(len=4)  :: cyear
   character(len=2)  :: cmonth
   character(len=40) :: fulldir,ftype
   type(hycomfile) :: hfile


   integer :: i, i2, i3,ilon,jlat,k ,k2
   real, pointer, dimension(:,:) :: csaln, cspres, ctemp, ctpres
   real, pointer, dimension(:,:) :: msaln , mtemp, msaln1, mtemp1, mpres1
   real, pointer, dimension(:,:) :: mpresd, msalnd, mtempd
   real, pointer, dimension(:,:) :: levsaln, levtemp

   ! Error flags - used locally by nersc
   logical, pointer, dimension(:,:) :: errorflag_sal,errorflag_temp

   ! invalid profile - per mersea definitions
   logical, pointer, dimension(:)   :: invalid_sal,invalid_temp
   real, pointer, dimension(:)   :: cslon,ctlon,cslat,ctlat,csjuld,ctjuld

   real :: lat_n,lon_n
   integer :: iprof, izlev, nprof, nzlev
   integer :: ipiv,jpiv
   integer :: validk
   integer :: iyear, imonth, iday
   logical :: l_ncdump=.false.
   integer :: fid_num, kdm
   character(len=10) :: cfid_num
   real, parameter :: undef999=-999

   
   ! Process input arguments
   if (iargc()/=2 .and. iargc()/=3 ) then
      print *,'Usage:'
      print *,'  argocmp argo-file hycom-file [ncdump]'
      stop
   else
      call getarg(1,corfile)
      call getarg(2,hyc_file)
      if (iargc()==3) then
         call getarg(3,tmpchar)
         if (trim(tmpchar)=='-ncdump') then
            l_ncdump=.true.
         else
            print *,'Usage:'
            print *,'  argocmp argo-file hycom-file [-ncdump]'
            stop
         end if
      end if
   end if

   ! Get bathymetry - lon/lat fields
   call xcspmd()
   call zaiost()
   call initconfmap(idm,jdm)
   call get_grid 

   ! What file is this? (daily? weekly? restart? pak?)
   ftype=getfiletype(trim(hyc_file))

   ! Initialize the hycomfile module - and retrieve vertical levels
   call initHF(hfile,trim(hyc_file),trim(ftype))
   kdm=vdim(hfile)

   ! Get the coriolis schtuff - for now we can deal with one input file only
   call read_coriolis(trim(corfile),csaln,cspres,cslon,cslat,csjuld,floatid, &
                      errorflag_sal, invalid_sal, 'PSAL',undef999)
   call read_coriolis(trim(corfile),ctemp,ctpres,ctlon,ctlat,ctjuld,floatid, &
                      errorflag_temp,invalid_temp,'TEMP',undef999)
   if (any(cspres/=ctpres)) then
      print *,'Temperature and salinity is defined at different levels. Fix me.'
      stop
   end if
   nzlev = size(csaln,1)
   nprof = size(csaln,2)


   ! Start working with hycom data file
   write (6,'(a)',advance='yes') 'Reading from '//trim(hyc_file)


   ! Start allocating fields now..
   allocate(dp       (idm,jdm))    ! Holds 2D vars on original grid
   allocate(pres     (idm,jdm))    ! Holds 2D vars on original grid
   allocate(presold  (idm,jdm))    ! Holds 2D vars on original grid
   allocate(temp     (idm,jdm))    ! Holds 2D vars on original grid
   allocate(saln     (idm,jdm))    ! Holds 2D vars on original grid
   allocate(msaln1(kdm,nprof))     ! Holds 2D vars on orig vertical grid pr station
   allocate(mtemp1(kdm,nprof))     ! Holds 2D vars on orig vertical grid pr station
   allocate(mpres1(kdm,nprof))     ! Holds 2D vars on orig vertical grid pr station
   allocate(msalnd(nzlev,nprof))   ! Holds 2D vars on argo vertical grid pr station
   allocate(mtempd(nzlev,nprof))   ! Holds 2D vars on argo vertical grid pr station
   allocate(mpresd(nzlev,nprof))   ! Holds 2D vars on argo vertical grid pr station
   allocate(msaln (nzlev,nprof))   ! Holds 2D vars on argo vertical grid pr station
   allocate(mtemp (nzlev,nprof))   ! Holds 2D vars on argo vertical grid pr station


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Start working through the original vertical model grid
   mpres1=0.
   pres=0.
   presold=0.
   allocate(valid_profile(nprof))
   valid_profile=.true.
   do k=1,kdm

      call HFReadDPField(hfile,dp,idm,jdm,k,1)
      call HFReadField(hfile,saln,idm,jdm,'saln    ',k,1)
      call HFReadField(hfile,temp,idm,jdm,'temp    ',k,1)
      pres=presold+dp/onem


      ! Put into model profile array 1. Same vertical levels as model
      ! Go through coriolis profiles -- find corresponding points in model
      do iprof=1,nprof

         ! Get closest point in model
         call oldtonew(cslat(iprof),cslon(iprof),lat_n,lon_n)
         call pivotp(lon_n,lat_n,ipiv,jpiv)

         ! Data point on grid
         if ( ipiv>1 .and. ipiv < idm .and. jpiv > 1 .and. jpiv < jdm .and. &
             depths(ipiv,jpiv)>1.) then
            if (k==1) then 
               print *,'model / obs point :',plon(ipiv,jpiv),plat(ipiv,jpiv), &
               cslon(iprof),cslat(iprof)
            end if
            mtemp1(k,iprof)=temp(ipiv,jpiv)
            msaln1(k,iprof)=saln(ipiv,jpiv)
            mpres1(k,iprof)=pres(ipiv,jpiv)
         else if (.not.(ipiv>1.and.ipiv<idm.and.jpiv>1.and.jpiv<jdm)) then
            if (k==1) print *,' obs point out of domain :', cslon(iprof),cslat(iprof)
            valid_profile(iprof)=.false.
         elseif (depths(ipiv,jpiv)<1.) then
            if (k==1) print *,' obs point on model land :', cslon(iprof),cslat(iprof)
            valid_profile(iprof)=.false.
         end if
      end do ! walk through stations
      presold=pres
   end do ! walk through levels


   ! We have the model stations on its own vertical levels. Interpolate to argo
   ! levels
   mtemp=undef
   msaln=undef
   do iprof=1,nprof

      ! 1) Get valid depth points from argo profile
      ! TODO - Use error flag
      validk=0
      do k=1,nzlev
         if (cspres(k,iprof)/=undef999) then
            validk=validk+1
            cspres(validk,iprof)=cspres(k,iprof)
            csaln (validk,iprof)=csaln (k,iprof)
            ctemp (validk,iprof)=ctemp (k,iprof)
         end if
      end do


      ! No valid points in the vertical
      if (validk==0) then
         valid_profile(iprof)=.false.

      ! Valid points in the vertical
      else if (valid_profile(iprof)) then

         if (validk<nzlev) then
            cspres(validk+1:nzlev,iprof)=undef999
         end if

         ! Initialize the spline calculations
         call spline_calc_ini_frominput('spline',cspres(1:validk,iprof),validk)

         ! Interpolate model to obs levels
         call spline_calc_1d(mtemp1(:,iprof),mpres1(:,iprof),mtemp(1:validk,iprof),validk,kdm)
         call spline_calc_1d(msaln1(:,iprof),mpres1(:,iprof),msaln(1:validk,iprof),validk,kdm)
         where (abs(mtemp-undef)<abs(undef)*1e-4) mtemp=undef999
         where (abs(msaln-undef)<abs(undef)*1e-4) msaln=undef999

         ! Release spline arrays for next profile (they may have different vertical levels)
         call clear_spline_calc()

         ! Also get discrete model layer values
         mtempd(:,iprof)=undef999
         msalnd(:,iprof)=undef999
         do k=1,validk
         do k2=1,kdm
            if (k2==1) then
               if ( mpres1(k2,iprof)>cspres(k,iprof)) then
                  mtempd(k,iprof)=mtemp1(k2,iprof)
                  msalnd(k,iprof)=msaln1(k2,iprof)
                  mpresd(k,iprof)=mpres1(k2,iprof)*.5
               end if
            else
               if ( mpres1(k2,iprof)>cspres(k,iprof) .and.  mpres1(k2-1,iprof)<cspres(k,iprof)) then
                  mtempd(k,iprof)=mtemp1(k2,iprof)
                  msalnd(k,iprof)=msaln1(k2,iprof)
                  mpresd(k,iprof)=(mpres1(k2,iprof) + mpres1(k2-1,iprof))*.5
               end if
            end if
         end do
         end do
      end if ! Valid points
   end do ! Loop  over stations


   ! Interpolation done - the rest is diagnostics


   ! Now copy the original file to a new one, then dump model values in the
   ! netcdf file
   if (l_ncdump) then

      call dump_coriolis(trim(corfile),'HYCOM_','PSAL',msaln,nzlev,nprof,undef999)
      call dump_coriolis(trim(corfile),'HYCOM_','TEMP',mtemp,nzlev,nprof,undef999)

      ! Read levitus data on its own levels 
      call levitus_setup('WOA2001')
      call levitus_interp('temperature',csjuld,cslon,cslat,cspres,levtemp,nprof,nzlev,undef999)
      call levitus_interp('salinity',   csjuld,cslon,cslat,cspres,levsaln,nprof,nzlev,undef999)
      !call stationsInterpLevitus('temperature',csjuld,cslon,cslat,cspres,levtemp,nprof,nzlev)
      !call stationsInterpLevitus('salinity'   ,csjuld,cslon,cslat,cspres,levsaln,nprof,nzlev)
      !modcorfile='TOPAZ_'//trim(corfile)
      !call system('cp '//corfile//' '//modcorfile)
      print *,'WOA2001 salinity '
      call dump_coriolis(trim(corfile),'WOA2001_','PSAL',levsaln,nzlev,nprof,undef999)
      print *,'WOA2001 temperature '
      call dump_coriolis(trim(corfile),'WOA2001_','TEMP',levtemp,nzlev,nprof,undef999)

      ! Read levitus data on its own levels 
      call levitus_setup('WOA2005')
      call levitus_interp('temperature',csjuld,cslon,cslat,cspres,levtemp,nprof,nzlev,undef999)
      call levitus_interp('salinity',   csjuld,cslon,cslat,cspres,levsaln,nprof,nzlev,undef999)
      !call stationsInterpLevitus('temperature',csjuld,cslon,cslat,cspres,levtemp,nprof,nzlev)
      !call stationsInterpLevitus('salinity'   ,csjuld,cslon,cslat,cspres,levsaln,nprof,nzlev)


      !modcorfile='TOPAZ_'//trim(corfile)
      !call system('cp '//corfile//' '//modcorfile)
      print *,'WOA2005 salinity '
      call dump_coriolis(trim(corfile),'WOA2005_','PSAL',levsaln,nzlev,nprof,undef999)
      print *,'WOA2005 temperature '
      call dump_coriolis(trim(corfile),'WOA2005_','TEMP',levtemp,nzlev,nprof,undef999)
   end if



   ! Now dump the profiles as text files - used by google map routines (and
   ! others)
   call system("mkdir -p ArgoCmp")
   where (errorflag_sal ) csaln=undef999
   where (errorflag_temp) ctemp=undef999
   do iprof=1,nprof
     if (valid_profile(iprof)) then

        write(cii,'(i3.3)') iprof

        !write(clon,'(f7.2)') cslon(iprof)     ; clon=adjustl(clon)
        !write(clat,'(f7.2)') cslat(iprof)     ; clat=adjustl(clat)
        write(clon,'(f7.2)') cslon(iprof)     ; clon=adjustl(clon)
        write(clat,'(f7.2)') cslat(iprof)     ; clat=adjustl(clat)

        write(cjuld,'(i10.10)') floor(csjuld(iprof))
        call juliantodate(floor(csjuld(iprof)) ,iyear,imonth,iday,1950,1,1)
        write(cdate,'(i4.4,i2.2,i2.2)') iyear,imonth,iday
        write(cyear,'(i4.4)') iyear
        write(cmonth,'(i2.2)') imonth

        ! Read float id number
        print '(a,a,3i5)','floatid , date is :',floatid(iprof)(:),iyear, imonth, iday
        read(floatid(iprof)(:),*,iostat=ios)   fid_num


        ! This neglects floats which are not all-digits....
        if (ios==0) then 
           write(cfid_num,'(i10.10)') fid_num

           ! KAL -- create new directory for files
           fulldir="ArgoCmp/"//cyear//"/"//cmonth
           call system("[ ! -d "//trim(fulldir)//" ] && mkdir -p "//trim(fulldir))

           !cprof='profile_date'//trim(cdate)//'_id'//trim(floatid(iprof)(:))//'_pos'//trim(clon)//'Ex'//trim(clat)//'N'
           cprof='profile_date'//trim(cdate)//'_id'//cfid_num//'_pos'//trim(clon)//'Ex'//trim(clat)//'N'
           !cprof='ArgoCmp/'//cprof
           cprof=trim(fulldir)//'/'//cprof


           ! KAL - TODO - use error flag here
           open(10,file=trim(cprof))
           do izlev=1,nzlev
           if (ctpres(izlev,iprof)/=undef999) then
              !print *,mtemp(izlev,iprof)
              write(10,'(f8.1,7f12.5)') ctpres(izlev,iprof),csaln(izlev,iprof),msaln(izlev,iprof), &
                                   ctemp(izlev,iprof),mtemp(izlev,iprof), &
                                   mtempd(izlev,iprof),msalnd(izlev,iprof),mpresd(izlev,iprof)
           end if
           end do
           close(10)

        else
           print *,'No ArgoCmp data produced for this profile -- argocmp does not handle non-digit profile ids (yet)'
        end if
     end if
  end do

   ! Dump a profile legend file as well
   open(10,file='ArgoCmp/profile_legend')
   write(10,*) 'Files contain the following data columnes for each profile'
   write(10,*) 'Column ' , ' Quantity'
   write(10,*) '1)     ' , ' Argo  pressure'
   write(10,*) '2)     ' , ' Argo  Salinity'
   write(10,*) '3)     ' , ' Model Salinity (in Argo point)'
   write(10,*) '4)     ' , ' Argo  Temperature'
   write(10,*) '5)     ' , ' Model Temperature (in Argo point)'
   write(10,*) '6)     ' , ' Model Temperature (in model layer center)'
   write(10,*) '7)     ' , ' Model Salinity    (in model layer center)'
   write(10,*) '8)     ' , ' Model Pressure    (in model layer center)'
   write(10,*)
   write(10,*) 'Filename Legend:'
   write(10,*) 'profile_date[date]_id[float id]_pos[position]'
   write(10,*) '[date] is  in format year month day of month (yyyymmdd)'
   write(10,*) '[float id] is float identifier (WMO float id)'
   write(10,*) 'position is in decimal degrees east x decimal degrees North'
   close(10)
   print *,'Profiles dumped  in directory "ArgoCmp"'




end program p_argocmp
