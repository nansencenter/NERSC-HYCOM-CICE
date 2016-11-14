module m_read_ncep

private :: humid, vapp

contains

   subroutine read_ncep(rt,plon,plat,depths)
      use mod_xc
      use mod_year_info 
      use mod_forcing_nersc
      use m_bilin_ncep_gauss
      !use m_bilin_ncep_gauss_mask
      use m_bilin_ecmwf2
      use m_ncvar_dims
      use m_ncvar_read
      use m_rotate
      use m_rotate_sphere
      use m_satvap
      use netcdf
      implicit none

      ! Field names
      character(len=*), parameter :: ncep_path='./NCEP/'

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

      type (year_info), intent(in) :: rt
      real, dimension(idm,jdm) :: plon,plat,depths
      character(len=5) :: cirec


      ! file names
      character(len=80) :: ficec,fprec,ftair,ftcdc,fuwnd, &
                           fvwnd, flmsk, fpres, fshum, fhgt

      integer , parameter :: lprt=6
      integer  :: i,j,k,ix,jy,irec
      logical, dimension(idm,jdm) :: modmsk, modmsk2

      real, dimension(:,:,:), allocatable :: &
         ncep_prec, ncep_pres, ncep_shum, ncep_tair, &
         ncep_tcdc, ncep_uwnd, ncep_vwnd, ncep_icec, &
         ncep_wspd, ncep_taux, ncep_tauy, ncep_rhum
      real, allocatable :: lon(:),lat(:),tmplat(:,:,:),tmplon(:,:,:),&
         tmpmsk(:,:,:) , ncep_lmsk(:,:), tmpv(:,:), ncep_hgt(:,:),   &
         tmphgt(:,:,:), lonr(:), latr(:) 

      real :: mlon(idm,jdm),tmp_icec(idm,jdm)
      real :: hgt(idm,jdm)
      real :: rfac,wndfac,cd_new,w4,wfact
      real :: dlon,dlat,dlonr,dlatr
      real :: vpair_x,vpair_w,vpair_i,humid_i,humid_w,sphum,vpair_1
      real :: svpair_w,svpair_i
      integer :: nlon,nlat,ndims,recdim,nlonr,nlatr
      integer , dimension(NF90_MAX_VAR_DIMS) :: dimsizes
      logical :: allocerr, ldump
      real :: satvapp, ax, bx,ncepcd,sathumid
      logical, save :: lfirst=.true.
      integer :: im1,ip1,jm1,jp1

      real, parameter :: & 
         aice     =9.5  ,&  ! --                 vapor pressure parameters
         bice     =7.66 ,&  ! k                    ..
         awater   =7.5  ,&  ! --                   ..
         bwater   =35.86    ! k                    ..


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
      ficec=ncep_path//'icec.sfc.gauss.'//rt%cyy//'.nc'
      fprec=ncep_path//'prate.sfc.gauss.'//rt%cyy//'.nc'
      ftair=ncep_path//'air.2m.gauss.'//rt%cyy//'.nc'
      ftcdc=ncep_path//'tcdc.eatm.gauss.'//rt%cyy//'.nc'
      fuwnd=ncep_path//'uwnd.10m.gauss.'//rt%cyy//'.nc'
      fvwnd=ncep_path//'vwnd.10m.gauss.'//rt%cyy//'.nc'
      flmsk=ncep_path//'land.sfc.gauss.nc'
      fhgt =ncep_path//'hgt.sfc.nc'
      fpres=ncep_path//'pres.sfc.gauss.'//rt%cyy//'.nc'
      fshum=ncep_path//'shum.2m.gauss.'//rt%cyy//'.nc'

      ! Check the allocation status of arrays
      allocerr=.false.
      if (lsynslp .and. .not.allocated(synslp)) then
         if (mnproc==1) write(lprt,*) 'SLP alloc error'
         allocerr=.true.
      else if (lsynrelhum .and. .not.allocated(synrelhum)) then
         if (mnproc==1) write(lprt,*) 'relhum alloc error' ;
         allocerr=.true.
      else if (lsynwnd .and. (         &
         .not.allocated(synuwind) .or. &
         .not.allocated(synvwind) .or. &
         .not.allocated(syntaux)  .or. &
         .not.allocated(syntauy)  .or. &
         .not.allocated(synwndspd) ) ) then
         if (mnproc==1) write(lprt,*) 'wind alloc error' ;
         allocerr=.true.
      else if (lsynairtmp .and. .not. allocated(synairtmp)) then
         if (mnproc==1) write(lprt,*) 'airtmp alloc error' ; 
         allocerr=.true.
      else if (lsynprecip .and. .not. allocated(synprecip)) then
         if (mnproc==1) write(lprt,*) 'precip alloc error' ; 
         allocerr=.true.
      end if
      if (allocerr) then
         if (mnproc==1) write(lprt,'(a)') 'Allocation error'
         call xcstop('(read_ncep)')
         stop '(read_ncep)'
      end if




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Get longitude  & latitude for data -- Two different grids 
      call ncvar_dims(flmsk,'lat',dimsizes,ndims,recdim)
      nlat=dimsizes(1) ; allocate(lat(nlat)) ; allocate(tmplat(nlat,1,1))
      call ncvar_read(flmsk,'lat', tmplat, nlat,1,1,1,1)
      lat=tmplat(:,1,1)
      !print *,'lat',lat
      deallocate(tmplat)

      call ncvar_dims(flmsk,'lon',dimsizes,ndims,recdim)
      nlon=dimsizes(1) ; allocate(lon(nlon)) ; allocate(tmplon(nlon,1,1))
      call ncvar_read(flmsk,'lon',tmplon,nlon,1,1,1,1)
      lon=tmplon(:,1,1)
      !print *,'lon',lon
      deallocate(tmplon)

      dlon=lon(2)-lon(1)
      do i=2,nlon
         if ( lon(i)-lon(i-1) /= dlon) then
            if (mnproc==1) write(lprt,'(a)') 'Longitude is not uniform !!'
            call xcstop('(forfun_ncep)')
            stop '(forfun_ncep)'
         end if
      end do



      ! Allocate temp fields
      allocate(ncep_icec(nlon,nlat,1))
      allocate(ncep_prec(nlon,nlat,1))
      allocate(ncep_tair(nlon,nlat,1))
      allocate(ncep_tcdc(nlon,nlat,1))
      allocate(ncep_uwnd(nlon,nlat,1))
      allocate(ncep_vwnd(nlon,nlat,1))
      allocate(ncep_wspd(nlon,nlat,1))
      allocate(ncep_taux(nlon,nlat,1))
      allocate(ncep_tauy(nlon,nlat,1))
      allocate(ncep_lmsk(nlon,nlat))
      allocate(ncep_shum(nlon,nlat,1))
      allocate(ncep_rhum(nlon,nlat,1))
      allocate(ncep_pres(nlon,nlat,1))
      allocate(tmpmsk(nlon,nlat,1))


      ! Depth mask
      modmsk=depths>.1



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NCEP land mask
      if (mnproc==1.and.ldump) write(lp, *) 'Using NCEP land mask'
      call ncvar_read(flmsk,vlmsk,tmpmsk   ,nlon,nlat, 1,1,1)
      ncep_lmsk=tmpmsk(:,:,1)

      !print *,'sum lmsk:',sum(ncep_lmsk(:,:))
      !do j=1,nlat
      !   print *,sum(ncep_lmsk(:,j))
      !end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NCEP geopot height
      call ncvar_dims(fhgt,'lat',dimsizes,ndims,recdim)
      nlatr=dimsizes(1) ; allocate(latr(nlatr)) ; allocate(tmplat(nlatr,1,1))
      call ncvar_read(fhgt,'lat', tmplat, nlatr,1,1,1,1)
      latr=tmplat(:,1,1)
      !print *,'lat',lat
      deallocate(tmplat)

      call ncvar_dims(fhgt,'lon',dimsizes,ndims,recdim)
      nlonr=dimsizes(1) ; allocate(lonr(nlonr)) ; allocate(tmplon(nlonr,1,1))
      call ncvar_read(fhgt,'lon',tmplon,nlonr,1,1,1,1)
      lonr=tmplon(:,1,1)
      !print *,'lon',lon
      deallocate(tmplon)
      allocate(ncep_hgt (nlonr,nlatr))

      allocate(tmphgt(nlonr,nlatr,1))
      if (mnproc==1.and.ldump) write(lp, *) 'Using NCEP geopot height (orography pressure correction'
      call ncvar_read(fhgt,vhgt,tmphgt   ,nlonr,nlatr, 1,1,1)
      ncep_hgt=tmphgt(:,:,1)
      dlonr=lonr(2)-lonr(1)
      dlatr=latr(2)-latr(1)

      ! Not ready yet ...
      !if (mnproc==1) write(lprt,'(a)') 'Fix bilinecmwf2r'
      !call xcstop('(read_ncep)')
      !stop '(read_ncep)'
      call bilin_ecmwf2(ncep_hgt,nlonr,nlatr,lonr(1),latr(1),dlonr,dlatr,&
         hgt,mlon,plat,depths)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NCEP ice concentration
      if (mnproc==1.and.ldump) write(lp, *) 'Using NCEP ice conc (involved in relhum)'
      call ncvar_read(ficec,vicec,ncep_icec,nlon,nlat,1,irec,irec)
      call bilin_ncep_gauss(ncep_icec(:,:,1),nlon,nlat,lon,lat,&
           tmp_icec(:,:),mlon,plat)





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NCEP precipitation   kg/(m^2 s) -> m/month (rho=1000)
      if (lsynprecip) then
         if (mnproc==1.and.ldump) write(lp, *) 'Using NCEP synoptic precipitation '
         call ncvar_read(fprec,vprec,ncep_prec,nlon,nlat,1,irec,irec)
         call bilin_ncep_gauss(ncep_prec(:,:,1),nlon,nlat,lon,lat,&
               synprecip,mlon,plat)

         ! convert prcipitation  from  kg/(m2*s) to m^3/m^2s
         rfac=1e-3
         synprecip=  rfac* synprecip
      end if




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NCEP air temperature   K -> K
      if (lsynairtmp) then
         !print *,trim(ftair),' ',rt%cyy,irec
         !print *,nlon,nlat
         if (mnproc==1.and.ldump) write(lp,'(a)')' Using NCEP synoptic temperature'
         call ncvar_read(ftair,vtair,ncep_tair,nlon,nlat,1,irec,irec)
         call bilin_ncep_gauss(ncep_tair(:,:,1),nlon,nlat,lon,lat,&
              synairtmp,mlon,plat)
      end if





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NCEP clouds   % -> []  (0-1)
      if (lsynairtmp) then
         if (mnproc==1.and.ldump) write(lp,'(a)') ' Using NCEP synoptic clouds '
         call ncvar_read(ftcdc,vtcdc,ncep_tcdc,nlon,nlat,1,irec,irec)
         call bilin_ncep_gauss(ncep_tcdc(:,:,1),nlon,nlat,lon,lat,&
              synclouds,mlon,plat)
         synclouds=0.01*synclouds
         synclouds=min(1.0,max(synclouds,0.))
         !print *,maxval(synclouds),minval(synclouds)
      end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NCEP winds
      if (lsynwnd) then
         if (mnproc==1.and.ldump) write(lp,'(a)') ' Using NCEP synoptic winds '

         ! Increase land mask to allow for interpolation below
         !modmsk2=.false.
         !marg_loc=nbdy-1
         !do j=1,jdm
         !do i=1,idm
         !   if (modmsk(i,j)) modmsk2(i-1:i+1,j-1:j+1) = modmsk(i,j)
         !end do
         !end do

         call ncvar_read(fuwnd,vuwnd,ncep_uwnd,nlon,nlat,1,irec,irec)
         call ncvar_read(fvwnd,vvwnd,ncep_vwnd,nlon,nlat,1,irec,irec)

         call bilin_ncep_gauss(ncep_uwnd(:,:,1),nlon,nlat,lon,lat,&
              synuwind,mlon,plat)
         call bilin_ncep_gauss(ncep_vwnd(:,:,1),nlon,nlat,lon,lat,&
              synvwind,mlon,plat)

         ! Rotate winds to model grid 
         call rotate(synuwind, synvwind,  plat    , mlon    ,idm,jdm,'l2m')
         !print *,'NCEP calls rotate_sphere'
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
         ! A word of 
         if (mnproc==1.and.ldump) write(lp,'(a)') ' Using NCEP synoptic SLP '
         call ncvar_read(fpres,vpres,ncep_pres,nlon,nlat,1,irec,irec)
         call bilin_ncep_gauss(ncep_pres(:,:,1),nlon,nlat,lon,lat,&
              synslp,mlon,plat)

         synslp(:,:)=0.01*synslp(:,:)    ! Convert from Pa to mbar 
         synslp(:,:)=synslp+hgt*9.8*1.3*0.01 ! Correction for geopot.
         if (ldump.and.mnproc==1) &
         print *,'max/min geopot hgt adj:',maxval(hgt)*9.8*1.3*0.01, &
           minval(hgt)*9.8*1.3*0.01




         if (mnproc==1.and.ldump) then
            write(lp,'(a)') ' Using NCEP synoptic specific humidity'
            call flush(lp)
         end if
         call ncvar_read(fshum,vshum,ncep_shum,nlon,nlat,1,irec,irec)
         call bilin_ncep_gauss(ncep_shum(:,:,1),nlon,nlat,lon,lat,&
              synrelhum,mlon,plat)
         if (mnproc==1.and.ldump) then
            write(lp,'(a)') ' End   NCEP synoptic specific humidity'
            call flush(lp)
         end if

         do jy=1,jdm
         do ix=1,idm

            ! NB -- no correction for air temp (due to slp corr above)
            !print *,i,j,tmp_icec(i,j)
    
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
            !print *,ix,jy,ncep_tair(ix,jy,1),sathumid,ncep_shum(ix,jy,1)
            synrelhum(ix,jy)=min(max(0.0,synrelhum(ix,jy)/sathumid),1.)

            !call xcstop('(check pressure etc in ncep fields)')
            !stop '(check pressure etc in ncep fields)'
         end do
         end do
      end if

      if (mnproc==1.and.ldump)  then
         print '(a)','####################################################'
         print *
      end if

      lfirst=.false.


      ! Clean up workspace
      deallocate(ncep_icec)
      deallocate(ncep_prec)
      deallocate(ncep_tair)
      deallocate(ncep_tcdc)
      deallocate(ncep_uwnd)
      deallocate(ncep_vwnd)
      deallocate(ncep_wspd)
      deallocate(ncep_taux)
      deallocate(ncep_tauy)
      deallocate(ncep_lmsk)
      deallocate(ncep_hgt)
      deallocate(tmpmsk)
      deallocate(tmphgt)
      deallocate(lat)
      deallocate(lon)
      deallocate(latr)
      deallocate(lonr)
      deallocate(ncep_rhum)
      deallocate(ncep_shum)
      deallocate(ncep_pres)



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
