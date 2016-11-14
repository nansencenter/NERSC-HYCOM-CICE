module m_read_ncep

private :: humid, vapp

contains

   subroutine read_ncep(rt,plon,plat,depths)
      use mod_xc
      use mod_year_info22
      use mod_forcing_nersc
      use m_bilin_ncep_gauss
      !use m_bilin_ncep_gauss_mask
      use m_bilin_ecmwf2
      use m_ncvar_dims
      use m_ncvar_read
      use m_rotate
      use netcdf
      implicit none

      ! Field names
      character(len=200), save :: ncep_path='./NCEP/'

      ! Variable names
      character(len=*), parameter :: vicec='icec'
      character(len=*), parameter :: vprec='prate'
      character(len=*), parameter :: vtair='air'
      character(len=*), parameter :: vtcdc='tcdc'
      character(len=*), parameter :: vuwnd='uwnd'
      character(len=*), parameter :: vvwnd='vwnd'
      character(len=*), parameter :: vlmsk='land'
      character(len=*), parameter :: vhgt ='hgt'
      character(len=*), parameter :: vpres='pres'
      character(len=*), parameter :: vshum='shum'
      character(len=*), parameter :: vssrd='dswrf'

      ! vapor pressure parameters
      real, parameter :: & 
         aice     =9.5  ,&  
         bice     =7.66 ,&  ! k                    ..
         awater   =7.5  ,&  ! --                   ..
         bwater   =35.86    ! k                    ..

      type (year_info), intent(in) :: rt
      real, dimension(idm,jdm) :: plon,plat,depths
      character(len=4) :: clyy

      ! file names
      character(len=200) :: cenv
      character(len=80) ::  &
         ficec,fprec,ftair,ftcdc,fuwnd, fvwnd, flmsk, fpres,  &
         fshum, fhgt, fssrd, fssrd2 

      real, dimension(:,:,:), allocatable :: &
         nceptmp, ncep_icec, ncep_lmsk, ncep_hgt, tmplon, tmplat
      real, allocatable :: lon(:),lat(:),tmpdswrf(:,:)

      real :: mlon(idm,jdm),tmp_icec(idm,jdm)
      real :: rfac,wndfac,cd_new,w4,wfact,dlon,dlat
      real :: satvapp, ax, bx,ncepcd,sathumid
      integer :: nlon,nlat,ndims,recdim 
      integer , dimension(NF90_MAX_VAR_DIMS) :: dimsizes
      integer :: im1,ip1,jm1,jp1, maxrec
      integer  :: i,j,k,ix,jy,irec
      logical :: allocerr, ldump
      logical, save :: lfirst=.true.
      integer*4, external :: system
      integer*4 :: ret


      ! Check NCEPR path on first pass
      if (lfirst) then
         call getenv('NCEP_PATH',cenv)
         if (trim(cenv)/='') then ! prefer this path if present
            ncep_path=trim(cenv)//'/'
         end if
         ret=system('[ -d '//trim(ncep_path)//' ]')
         if (ret /=0 ) then
            print *
            print *,'The directory  '//trim(ncep_path)//' does not exist.'
            print *,'Make sure that '
            print *,' - You have linked the NCEP data into this catalogue, and '
            print *,'   that the variable NCEP_PATH is empty'
            print *,' - OR, the variable NCEP_PATH is set to the location of '
            print *,'   the ERA40 data'
            call exit(1)
         end if
      end if
      !call xcstop('read_ncep not changed to using tiles')
      !stop 'read_ncep not changed to using tiles'
      mlon = plon ; where (mlon < 0.) mlon=mlon+360.

      ! We are reading record number :
      irec = rt%ihh/6 + rt%idd*4 +1

      ldump=lfirst.or.mod(irec,40)==0

      if (mnproc==1.and.ldump) then
         print *
         write(lp,'(a)')'####################################################'
         write(lp,'(a)')'Reading NCEP synoptic forcing fields '
         write(lp,'(a,i5,a,i4,a,i3,a,i2)')'Record=',irec,' Year=',rt%iyy, &
            ' Day=',rt%idd,' hour=',rt%ihh
      else
         write(lp,'(a,i5,a,i4,a,i3,a,i2,a,i2,a,i3)')'Record=',irec,' Year=',rt%iyy, &
            ' Day=',rt%idd,' hour=',rt%ihh, ' month=',rt%imm,' dayinmonth=',rt%idm+1
      end if



      ! Which file to open
      write(clyy,'(i4.4)') rt%iyy-1
      ficec=trim(ncep_path)//'icec.sfc.gauss.'//rt%cyy//'.nc'
      fprec=trim(ncep_path)//'prate.sfc.gauss.'//rt%cyy//'.nc'
      ftair=trim(ncep_path)//'air.2m.gauss.'//rt%cyy//'.nc'
      ftcdc=trim(ncep_path)//'tcdc.eatm.gauss.'//rt%cyy//'.nc'
      fuwnd=trim(ncep_path)//'uwnd.10m.gauss.'//rt%cyy//'.nc'
      fvwnd=trim(ncep_path)//'vwnd.10m.gauss.'//rt%cyy//'.nc'
      fpres=trim(ncep_path)//'pres.sfc.gauss.'//rt%cyy//'.nc'
      fshum=trim(ncep_path)//'shum.2m.gauss.'//rt%cyy//'.nc'
      fssrd=trim(ncep_path)//'dswrf.sfc.gauss.'//rt%cyy//'.nc'
      fssrd2=trim(ncep_path)//'dswrf.sfc.gauss.'//clyy//'.nc'

      flmsk=trim(ncep_path)//'land.sfc.gauss.nc'
      !fhgt =trim(ncep_path)//'hgt.sfc.nc' ! Old  on regular grid
      fhgt =trim(ncep_path)//'hgt.sfc.gauss.nc' ! New - on gauss grid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Get longitude  & latitude for data -- Two different grids 
      call ncvar_dims(flmsk,'lat',dimsizes,ndims,recdim)
      nlat=dimsizes(1) ; allocate(lat(nlat)) ; allocate(tmplat(nlat,1,1))
      call ncvar_read(flmsk,'lat', tmplat, nlat,1,1,1,1)
      lat=tmplat(:,1,1)
      deallocate(tmplat)

      call ncvar_dims(flmsk,'lon',dimsizes,ndims,recdim)
      nlon=dimsizes(1) ; allocate(lon(nlon)) ; allocate(tmplon(nlon,1,1))
      call ncvar_read(flmsk,'lon',tmplon,nlon,1,1,1,1)
      lon=tmplon(:,1,1)
      deallocate(tmplon)

      dlon=lon(2)-lon(1)
      do i=2,nlon
         if ( lon(i)-lon(i-1) /= dlon) then
            if (mnproc==1) write(lp,'(a)') 'Longitude is not uniform !!'
            call xcstop('(forfun_ncep)')
         end if
      end do

      ! Allocate temp fields
      allocate(ncep_icec (nlon,nlat,1))
      allocate(ncep_lmsk (nlon,nlat,1))
      allocate(ncep_hgt  (nlon,nlat,1))
      allocate(nceptmp   (nlon,nlat,1))
      allocate(tmpdswrf  (nlon,nlat))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NCEP land mask
      if (mnproc==1.and.ldump) write(lp, *) 'Using NCEP land mask'
      call ncvar_read(flmsk,vlmsk,ncep_lmsk   ,nlon,nlat, 1,1,1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NCEP geopot height
      call ncvar_read(fhgt,vhgt,ncep_hgt   ,nlon,nlat, 1,1,1)
      if (mnproc==1.and.ldump) write(lp, *) 'Using NCEP geopot height (orography pressure correction'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NCEP ice concentration
      if (mnproc==1.and.ldump) write(lp, *) 'Using NCEP ice conc (involved in relhum)'
      call ncvar_read(ficec,vicec,ncep_icec,nlon,nlat,1,irec,irec)
      call bilin_ncep_gauss(ncep_icec(:,:,1),nlon,nlat,lon,lat, tmp_icec(:,:),mlon,plat)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NCEP precipitation   kg/(m^2 s) -> m/month (rho=1000)
      if (lsynprecip) then
         if (mnproc==1.and.ldump) write(lp, *) 'Using NCEP synoptic precipitation '
         call ncvar_read(fprec,vprec,nceptmp,nlon,nlat,1,irec,irec)
         call bilin_ncep_gauss(nceptmp(:,:,1),nlon,nlat,lon,lat,&
               synprecip,mlon,plat)

         ! convert prcipitation  from  kg/(m2*s) to m^3/m^2s
         synprecip=  synprecip*1e-3
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NCEP air temperature   K -> K
      if (lsynairtmp) then
         if (mnproc==1.and.ldump) write(lp,'(a)')' Using NCEP synoptic temperature'
         call ncvar_read(ftair,vtair,nceptmp,nlon,nlat,1,irec,irec)
         call bilin_ncep_gauss(nceptmp(:,:,1),nlon,nlat,lon,lat,&
              synairtmp,mlon,plat)
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NCEP clouds   % -> []  (0-1)
      if (lsynairtmp) then
         if (mnproc==1.and.ldump) write(lp,'(a)') ' Using NCEP synoptic clouds '
         call ncvar_read(ftcdc,vtcdc,nceptmp,nlon,nlat,1,irec,irec)
         call bilin_ncep_gauss(nceptmp(:,:,1),nlon,nlat,lon,lat,&
              synclouds,mlon,plat)
         synclouds=0.01*synclouds
         synclouds=min(1.0,max(synclouds,0.))
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NCEP winds
      if (lsynwnd) then
         if (mnproc==1.and.ldump) write(lp,'(a)') ' Using NCEP synoptic winds '
         call ncvar_read(fuwnd,vuwnd,nceptmp,nlon,nlat,1,irec,irec)
         call bilin_ncep_gauss(nceptmp(:,:,1),nlon,nlat,lon,lat,&
              synuwind,mlon,plat)
         call ncvar_read(fvwnd,vvwnd,nceptmp,nlon,nlat,1,irec,irec)
         call bilin_ncep_gauss(nceptmp(:,:,1),nlon,nlat,lon,lat,&
              synvwind,mlon,plat)

         ! Rotate winds to model grid 
         call rotate(synuwind, synvwind,  plat    , mlon    ,idm,jdm,'l2m')
         !call rotate_sphere(synuwind, synvwind,  plat    , mlon    ,idm,jdm,'l2m')

         ! Derived  wind stress and absolute windspeed
         ncepcd=0.0012

         !$OMP PARALLEL DO PRIVATE(j,i,wndfac,w4,wfact,cd_new,im1,jm1,ip1,jp1) &
         !$OMP SCHEDULE(STATIC,jblk)
         do j=1,jdm
         do i=1,idm
            im1=max(1,i-1)
            ip1=min(idm,i+1)
            jm1=max(1,j-1)
            jp1=min(jdm,j+1)

            synwndspd(i,j)=sqrt(synuwind(i,j)*synuwind(i,j)+ &
                               synvwind(i,j)*synvwind(i,j))

            wndfac=(1.+sign(1.,synwndspd(i,j)-11.))*.5
            cd_new=(0.49+0.065*synwndspd(i,j))*1.0e-3*wndfac+cd*(1.-wndfac)

            w4   =.25*(synvwind(im1,jp1)+synvwind(i,jp1)+ &
                       synvwind(im1,j  )+synvwind(i,j  ))
            wfact=sqrt(synuwind(i,j)*synuwind(i,j)+w4*w4)*airdns*cd_new
            syntaux(i,j)=synuwind(i,j)*wfact

            w4   =.25*(synuwind(i  ,jm1)+synuwind(ip1,jm1)+&
                       synuwind(ip1,j  )+synuwind(i  ,j  ))
            wfact=sqrt(synvwind(i,j)*synvwind(i,j)+w4*w4)*airdns*cd_new
            syntauy(i,j)=synvwind(i,j)*wfact
         end do
         end do
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!NCEP relative humidity and pressure
      if (lsynslp .or. lsynrelhum) then
         if (mnproc==1.and.ldump) write(lp,'(a)') ' Using NCEP synoptic SLP '
         call ncvar_read(fpres,vpres,nceptmp,nlon,nlat,1,irec,irec)
         nceptmp=0.01*nceptmp                  ! Convert from Pa to mbar 
         nceptmp=nceptmp+ncep_hgt*9.8*1.3*0.01 ! Simple correction for geopot. height
         call bilin_ncep_gauss(nceptmp(:,:,1),nlon,nlat,lon,lat,&
              synslp,mlon,plat)

         if (ldump.and.mnproc==1) &
         print *,'max/min geopot hgt adj:',maxval(ncep_hgt)*9.8*1.3*0.01, &
           minval(ncep_hgt)*9.8*1.3*0.01

         if (mnproc==1.and.ldump) then
            write(lp,'(a)') ' Using NCEP synoptic specific humidity'
            call flush(lp)
         end if
         call ncvar_read(fshum,vshum,nceptmp,nlon,nlat,1,irec,irec)
         call bilin_ncep_gauss(nceptmp(:,:,1),nlon,nlat,lon,lat,&
              synrelhum,mlon,plat)
         do jy=1,jdm
         do ix=1,idm
            if (tmp_icec(ix,jy)>.2) then
               ax=aice
               bx=bice
            else
               ax=awater
               bx=bwater
            end if

            ! Saturation Vapour pressure (sea level)
            satvapp=vapp(ax,bx,synairtmp(ix,jy))

            ! Saturation humidity
            sathumid=humid(synslp(ix,jy)*100.,satvapp)

            ! Relative humidty
            synrelhum(ix,jy)=min(max(0.0,synrelhum(ix,jy)/sathumid),1.)
         end do
         end do
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!NCEP downwelling shortwave radiation
      if (lsynshwflx ) then
         if (mnproc==1.and.ldump) write(lp,'(a)') ' Using NCEP synoptic DSWRF '
         print *,'NB: DSWRF not testet thoroughly'
         call ncvar_read(fssrd,vssrd,nceptmp,nlon,nlat,1,irec,irec)
         tmpdswrf(:,:)=nceptmp(:,:,1)

         ! ssrd is averaged over 6 next hrs, read last ssrd field as well
         if (irec==1)  then
            call ncvar_dims(fssrd2,vssrd,dimsizes,ndims,recdim)
            maxrec=dimsizes(recdim)
            call ncvar_read(fssrd2,vssrd,nceptmp,nlon,nlat,1,maxrec,maxrec)
         else 
            call ncvar_read(fssrd ,vssrd,nceptmp,nlon,nlat,1,irec-1,irec-1)
         end if
         nceptmp(:,:,1)=nceptmp(:,:,1)+tmpdswrf
         nceptmp=nceptmp*.5
         call bilin_ncep_gauss(nceptmp(:,:,1),nlon,nlat,lon,lat,&
              synshwflx,mlon,plat)
      end if


      if (mnproc==1.and.ldump)  then
         print '(a)','####################################################'
         print *
      end if

      lfirst=.false.


      ! Clean up workspace
      deallocate(ncep_icec)
      deallocate(ncep_lmsk)
      deallocate(ncep_hgt)
      deallocate(nceptmp)
      deallocate(tmpdswrf)
      deallocate(lat)
      deallocate(lon)
   END subroutine read_ncep

   real function vapp(ax,bx,tx)
   implicit none
     real, intent(in) :: ax,bx,tx
      vapp=611.*10.**(ax*(tx-273.16)/(tx-bx))
   end function vapp

   real function humid(px,ex)
   implicit none
     real, intent(in) :: px,ex
      humid=.622*ex/(px-.378*ex)
   end function humid

end module m_read_ncep
