module mod_storm

contains

   subroutine read_storm(rt,rtinit,plon,plat,depths,bailout)
      use mod_xc
      use mod_year_info 
      use mod_forcing_nersc
      !use mod_common_ice
      use m_ncvar_dims
      use m_ncvar_read
      use m_rotate
      use m_relhumid
      use m_satvap
      use netcdf
      implicit none

      ! Data catalogue name
      !character(len=*), parameter :: era40_path='./ERA40-1.125x1.125/'
      character(len=*), parameter :: storm_path='./Storm/'

      ! Variable names - Storm version (not really needed)

      type (year_info), intent(in) :: rt,rtinit
      real, dimension(idm,jdm), intent(in) :: plon,plat,depths
      logical, intent(inout) :: bailout
      character(len=5) :: cirec


      ! File names
      character(len=80) :: fprec, ftair, ftcc, fuwnd, fvwnd, fblh, &
         flmsk, fpres, fdewp, fssrd, fssrd2, ficec,fz
      character(len=25) :: cdate,cpost
      character(len= 2) :: cgap

      integer  :: i,j,k,ix,jy,irec,maxrec,margin_loc,irec_dummy

      real, dimension(:,:,:), allocatable :: &
         storm_prec, storm_pres, storm_dewp, storm_tair, &
         storm_tcc,  storm_uwnd, storm_vwnd, storm_ssrd
      real, allocatable ::  storm_lmsk(:,:,:), storm_rhum(:,:), tmpfld(:,:), rdfld(:,:), &
         dglon(:,:),dglat(:,:)

      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: mlon
      real :: rfac,wndfac,cd_new,w4,wfact,svpair, svpdew
      real :: dlon,dlat,flon,flat,llon,llat,dummy(1,1,1)
      real :: vpair_x,vpair_w,vpair_i,humid_i,humid_w,sphum,vpair_1
      integer :: sidm,sjdm,ndims,recdim
      integer , dimension(NF90_MAX_VAR_DIMS) :: dimsizes
      real :: satvapp, ax, bx,era40cd,sathumid
      integer :: ndx,ndy
      logical :: lfirst=.true.
      integer :: hgap
      integer :: im1,jm1,ip1,jp1

      character(len=4)  :: cyy

      real, parameter :: t0deg  =273.15



      !if (mnproc==1) print *,'ERA40 -- add time stamp check!'

      !call xcstop('read_rea40 not changed to using tiles')
      !stop 'read_rea40 not changed to using tiles'

      ! We are reading record number [era40 starts on 1.Jan 00:00]
      ! Should have a more robust checking here ...
      irec_dummy  = rt%ihh/6 + rt%idd*4 +1
      irec=1

      ! Diagnostic output
      if (mod(irec_dummy,10)==0.or.lfirst) then
         if (mnproc==1) then
            WRITE(lp,*)
            if (lfirst) then
               write(lp,'(a,i6,a,i4,x,i3,x,i2)',advance='yes')'STORM/MM5 data from record: ',irec,&
               '  rt= ',rt%iyy,rt%idd,rt%ihh
            else
               write(lp,'(a,i6,a,i4,x,i3,x,i2)',advance='no') 'STORM/MM5 data from record: ',irec,&
               '  rt= ',rt%iyy,rt%idd,rt%ihh
            end if
         end if
     else
        if (mnproc==1) write(lp,'(a1)',advance='no')'.'
      end if

      mlon = plon ; where (mlon < 0.) mlon=mlon+360.

      ! Set file names
      write(cyy,'(i4.4)') rt%iyy+1
      ! Templated
      ! fuwnd=storm_path//'U10'//rt%cyy//rt%cmm//rt%cdm//'ACCUM_TIME'//'.grd.
      !fuwnd=storm_path//'rtot6h:'//rt%cyy//rt%cdm//'.10U.nc'
      !fvwnd=storm_path//'ans.6h.'//rt%cyy//'.10V.nc'
      !fdewp=storm_path//'ans.6h.'//rt%cyy//'.2D.nc'
      !ftair=storm_path//'ans.6h.'//rt%cyy//'.2T.nc'
      !fpres=storm_path//'ans.6h.'//rt%cyy//'.MSL.nc'
      !fssrd=storm_path//'fcs.6h.'//rt%cyy//'.SSRD.nc' ! Forecast field !
      !ftcc =storm_path//'ans.6h.'//rt%cyy//'.TCC.nc'


      !Calculate days between initial and this date
      hgap=datetojulian(rt%iyy,rt%imm,rt%idm,rtinit%iyy,rtinit%imm,rtinit%idm)

      ! Gap in hours
      hgap=hgap*24
      hgap=hgap-rtinit%ihh
      hgap=hgap+rt%ihh

      ! Encode gap in text
      !print *,hgap

      ! For now, restrict gap to <=72 hours
      if (hgap>72) then
         if (mnproc==1) then 
            print *,'Too big time gap for Storm forcing, bailing out'
         end if
         !call xcstop('(mod_storm)')
         !stop '(mod_storm)'
         bailout=.true.
         return
      else if (hgap==0) then
         !something fishy with first field...
         hgap=1
      end if

      write(cgap,'(i2.2)') hgap


      ! íntial date in text
      cdate=rtinit%cyy//rtinit%cmm//rtinit%cdm//rtinit%chh
      write(cdate(7:8),'(i2.2)') rtinit%idm+1
      cdate=trim(cdate)//'_'//cgap
      cpost='_'//trim(cdate)//'.grd'

      


      ! File names
      fprec=storm_path//'rtot6h'//trim(cpost)
      ftair=storm_path//'T2'    //trim(cpost)
      ftcc =storm_path//'ccov'  //trim(cpost)
      fuwnd=storm_path//'U10'   //trim(cpost)
      fvwnd=storm_path//'V10'   //trim(cpost)
      fdewp=storm_path//'sdpc'  //trim(cpost)
      fpres=storm_path//'slp'   //trim(cpost)
      fssrd=storm_path//'swdown'//trim(cpost)


      ! Check the allocation status of arrays
      if ((.not.allocated(synslp   ).and. lsynslp   ) .or. &
          (.not.allocated(synrelhum).and. lsynrelhum) .or. &
          (.not.allocated(synuwind ).and. lsynwnd   ) .or. &
          (.not.allocated(synvwind ).and. lsynwnd   ) .or. &
          (.not.allocated(syntaux  ).and. lsynwnd   ) .or. &
          (.not.allocated(syntauy  ).and. lsynwnd   ) .or. &
          (.not.allocated(synwndspd).and. lsynwnd   ) .or. &
          (.not.allocated(synairtmp).and. lsynairtmp) .or. &
          (.not.allocated(synprecip).and. lsynprecip) .or. &
          (.not.allocated(synshwflx).and. lsynshwflx)) then
         if (mnproc==1) write(lp,'(a)') 'ERA40 Allocation error'
         call xcstop('(read_storm)')
         stop '(read_storm)'
      end if




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Get longitude  & latitude for data  - use land mask file
      call ncvar_dims(ftair,'z',dimsizes,ndims,recdim) ; 
      !print *,dimsizes(1:ndims)
      !sidm=dimsizes(2)
      !sjdm=dimsizes(1)
      sjdm=dimsizes(2)
      sidm=dimsizes(1)


      ! Allocate temp fields
      allocate(storm_prec(sidm,sjdm,1))
      allocate(storm_tair(sidm,sjdm,1))
      allocate(storm_tcc (sidm,sjdm,1))
      allocate(storm_uwnd(sidm,sjdm,1))
      allocate(storm_vwnd(sidm,sjdm,1))
      allocate(storm_pres(sidm,sjdm,1))
      allocate(storm_ssrd(sidm,sjdm,1))
      allocate(storm_dewp(sidm,sjdm,1))
      allocate(storm_lmsk(sidm,sjdm,1))
      allocate(storm_rhum(sidm,sjdm))
      allocate(tmpfld    (sidm,sjdm))
      allocate(dglon     (sidm,sjdm))
      allocate(dglat     (sidm,sjdm))
      !allocate(rdfld     (sjdm,sidm))



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! STORM precipitation   m/ 6hr -> m/s (rho=1000)
      if (lsynprecip) then
         if (mnproc==1.and.lfirst) write(lp,'(a)') ' Using STORM precipitation '
         call ncvar_read(fprec,'z',storm_prec,sidm,sjdm,1,irec,irec)
         call bilin_storm(storm_prec,sidm,sjdm,'test', synprecip,mlon,plat,depths)   

         ! Precipitation is accumulated field... mm over 6 hours
         synprecip=synprecip*1e-3/(6.*3600.)

         !print *,'precip: ',maxval(synprecip),minval(synprecip)

      end if




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! STORM air temperature   K 
      if (lsynairtmp) then
         if (mnproc==1.and.lfirst) write(lp,'(a)') ' Using STORM temperature'
         call ncvar_read(ftair,'z',storm_tair,sidm,sjdm,1,irec,irec)
         call bilin_storm(storm_tair,sidm,sjdm,'test', synairtmp,mlon,plat,depths)   
         !print *,minval(storm_tair),maxval(storm_tair)
         !print *,minval(synairtmp),maxval(synairtmp)
         !stop
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! STORM clouds   []  (0-1)
     if (lsynclouds) then
        if (mnproc==1.and.lfirst) write(lp,'(a)')' Using STORM clouds '
        call ncvar_read(ftcc,'z',storm_tcc,sidm,sjdm,1,irec,irec)
        call bilin_storm(storm_tcc ,sidm,sjdm,'test',synclouds,mlon,plat,depths)   
        
        synclouds=synclouds/100. ! % -> fractions [0-1]
     end if

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! STORM winds
      if (lsynwnd) then
         if (mnproc==1.and.lfirst) write(lp,'(a)')' Using STORM winds '

         call ncvar_read(fuwnd,'z',storm_uwnd,sidm,sjdm,1,irec,irec)
         call ncvar_read(fvwnd,'z',storm_vwnd,sidm,sjdm,1,irec,irec)

         ! Get lon/lat for data grid
         do j=1,sjdm
         do i=1,sidm
            call storm_xy2ll(real(i),real(j),'test',dglon(i,j),dglat(i,j))
         end do
         end do

         ! Rotate velocities on data grid to velocities on east/north direction
         !print *,maxval(abs(storm_uwnd))
         call rotate(storm_uwnd, storm_vwnd, dglat     , dglon     ,sidm,sjdm,'m2l')
         !print *,maxval(abs(storm_uwnd))


         ! Bilinear interpolation of new fields
         call bilin_storm(storm_uwnd ,sidm, sjdm, 'test', synuwind,mlon,plat,depths)   
         call bilin_storm(storm_vwnd ,sidm, sjdm, 'test', synvwind,mlon,plat,depths)   

         !if (mnproc==1) then
         !open(10,file='teststormuv')
         !do j=1,jdm
         !do i=1,idm
         !   write(10,*) plon(i,j),plat(i,j),synuwind(i,j),synvwind(i,j)
         !end do
         !end do
         !close(10)
         !end if


         ! Rotate winds to model grid 
         !print *,'Fix winds for STORM (rotation)'
         call rotate(synuwind,synvwind,plat,mlon,idm,jdm,'l2m')

         !if (mnproc==1) then
         !open(10,file='teststormuv2')
         !do j=1,jdm
         !do i=1,idm
         !   write(10,*) i,j,synuwind(i,j),synvwind(i,j),depths(i,j)
         !end do
         !end do
         !close(10)
         !end if
         !stop

         ! Derived  wind stress and absolute windspeed
         era40cd=0.0012

         !$OMP PARALLEL DO PRIVATE(j,i,wndfac,w4,wfact,cd_new) &
         !$OMP SCHEDULE(STATIC,jblk)
         do j=1,jdm
         do i=1,idm

            ip1=max(1,min(idm,i+1))
            jp1=max(1,min(jdm,j+1))
            im1=max(1,min(idm,im1))
            jm1=max(1,min(jdm,jm1))

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
! Pressure

         
      if (lsynslp) then
         if (mnproc==1.and.lfirst) write(lp,'(a)') ' Using STORM SLP '
         call ncvar_read(fpres,'z',storm_pres,sidm,sjdm,1,irec,irec)
         storm_pres=storm_pres*100. ! HPa -> Pa (needed later)
         call bilin_storm(storm_pres ,sidm, sjdm, 'test', synslp,mlon,plat,depths)   
         synslp=synslp*0.01 ! Pa -> mBar (HPa)
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Relative humidity
      if (lsynrelhum) then

         if (mnproc==1.and.lfirst) write(lp,'(a)') ' Using ERA40-derived relative humidity'
         call ncvar_read(fdewp,'z',storm_dewp,sidm,sjdm,1,irec,irec)

         ! SLP is needed to calc relative humidity (alternative - use slp0)
         if (.not.lsynslp) then
            write(lp,'(a)') ' ERROR: STORM slp is needed to calculate relhum ...  '
            call xcstop ('(read_storm)')
            stop '(read_storm)'
         end if

         ! relative humidity
         !$OMP PARALLEL DO PRIVATE (i,j,svpair,svpdew)
         do j=1,sjdm
         do i=1,sidm
            svpair=satvap(storm_tair(i,j,1))
            svpdew=satvap(storm_dewp(i,j,1)+t0deg) ! dewp in C
            storm_rhum(i,j)=relhumid(svpair,svpdew,storm_pres(i,j,1))*0.01

            ! Ahem ...
            storm_rhum(i,j)  =max(0.0,min(storm_rhum(i,j),1.0))

         enddo
         enddo
         !$OMP END PARALLEL DO
         !print *,maxval(storm_rhum)
         !print *,maxval(storm_tair)
         !print *,maxval(storm_dewp)+t0deg
         call bilin_storm(storm_rhum,sidm, sjdm, 'test', synrelhum,mlon,plat,depths)   
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SW radiation

      if (lsynssr) then
         if (mnproc==1.and.lfirst) write(lp,'(a)') ' Using STORM shortwave radiation '

         !call ncvar_read(fssrd,vssrd,era40_ssrd,nlon,nlat,1,irec,irec)
         call ncvar_read(fssrd,'z',storm_ssrd,sidm,sjdm,1,irec,irec)
         tmpfld=storm_ssrd(:,:,1)

!         ! ssrd is accumulated over 6hrs, midpoint is 6hr-3hr
!         if (irec==maxrec)  then
!            call ncvar_read(fssrd2,vssrd,era40_ssrd,nlon,nlat,1,1,1)
!         else if (irec<maxrec) then
!            call ncvar_read(fssrd ,vssrd,era40_ssrd,nlon,nlat,1,irec+1,irec+1)
!         else
!            print *,'You shouldnt be here ...'
!            print *,irec,maxrec
!            call xcstop('(read_era40)')
!            stop '(read_era40)'
!         end if
!         era40_ssrd(:,:,1)=era40_ssrd(:,:,1)+tmpfld
!         era40_ssrd=era40_ssrd/2.
!
!
         call bilin_storm(storm_ssrd ,sidm, sjdm, 'test', synshwflx,mlon,plat,depths)   
         !synshwflx=synshwflx/(3*3600) ! Wm^-2 * s - > Wm^-2
      end if


      ! Clean up workspace
      deallocate(storm_prec)
      deallocate(storm_tair)
      deallocate(storm_tcc )
      deallocate(storm_uwnd)
      deallocate(storm_vwnd)
      deallocate(storm_pres)
      deallocate(storm_ssrd)
      deallocate(storm_dewp)
      deallocate(storm_lmsk)
      deallocate(storm_rhum)
      deallocate(tmpfld)

      lfirst=.false.


   END subroutine read_storm


   subroutine bilin_storm(old,onx,ony,grid_id,new,newlon,newlat,newmask)   
   use mod_xc
   implicit none
   integer, intent(in) :: onx         ! x-dimension of old field
   integer, intent(in) :: ony         ! y-dimension of old field
   real, intent(in)    :: old(onx,ony)! old grid
   real, intent(out)   :: new   (idm,jdm) ! Interpolated field
   real, intent(in)    :: newlon(idm,jdm) ! Longitudes for new grid
   real, intent(in)    :: newlat(idm,jdm) ! Latitudes for new grid
   real, intent(in)    :: newmask(idm,jdm) ! Latitudes for new grid
   character(len=*), intent(in) ::  grid_id


   integer i,j,ib,jb
   integer ipos    !  index of i-pivot grid point in old grid
   integer jpos    !  index of j-pivot grid point in old grid
   real :: ripos,rjpos

   real aa,bb,a1,a2,a3,a4,lon
   integer :: numerr,sumerr
   character(len=2) tag2


   sumerr=0
   !open(10,file='blah',status='replace')
   ! Start interpolation
   !$OMP PARALLEL DO PRIVATE(i,j,ipos,jpos,lon,ib,jb,aa,bb,a1,a2,a3,a4) &
   !$OMP SCHEDULE(STATIC,jblk) REDUCTION(+:numerr)
   do j=1,jdm
   do i=1,idm
      numerr=0

      call storm_ll2xy(newlon(i,j),newlat(i,j),'test',ripos,rjpos)
      ipos=floor(ripos)
      jpos=floor(rjpos)

      !write(10,*) ipos,jpos,newlon(i,j),newlat(i,j)



      ! Error in case out of data domain
      numerr = numerr + (1-sign(1,jpos-1))/2
      numerr = numerr + (1-sign(1,ipos-1))/2
      numerr = numerr + (1-sign(1,onx-(ipos+1)))/2
      numerr = numerr + (1-sign(1,ony-(jpos+1)))/2
      !stop




      aa=ripos-ipos
      bb=rjpos-jpos
      ib=ipos+1
      jb=jpos+1
      !print *,numerr,ib,jb,onx,ony


      ! Error in case of wrong aa or bb
      numerr = numerr + abs(floor(aa))
      numerr = numerr + abs(floor(bb))

      sumerr=sumerr+numerr

      a1=(1.0-aa)*(1.0-bb)
      a2=aa*(1.0-bb)
      a3=aa*bb
      a4=(1.0-aa)*bb

      if (numerr==0) then
         new(i,j) = a1*old(ipos,jpos)+a2*old(ib,jpos)+a3*old(ib,jb)+a4*old(ipos,jb)
      else if (numerr/=0) then
         ! 1e20 is the "mask" value used in hycom
         if (newmask(i,j)<.1 .or. newmask(i,j)>1e20) then
            new(i,j) = 0.
         else
            new(i,j) = 0.
            print *,i,j,ipos,jpos,newlon(i,j),newlat(i,j)
            print*, '(Encountered ocean fields outside of forcing grid !)'
            stop '(bilin_storm)'
         end if
      endif
      !new(i,j)=jpos

      !if (ipos>1 .and.jpos>1.and.ipos<onx.and.jpos<ony) print *,ipos,jpos,numerr,aa,bb

   enddo
   enddo
   !$OMP END PARALLEL DO

end subroutine 

   subroutine storm_ll2xy(lon,lat,grid_id,x,y)
      implicit none
      real, parameter :: dskmc=12.  !Gitteravstand i km
      real, parameter :: ydim=230.  !y dimensions of mother domain
      real, parameter :: xdim=260.  !x dimensions of mother domain
      real, parameter :: true1=60.  ! true latitude of polar stereographic grid --- standardized for mm5
      real, parameter :: true2=60.
      real, parameter :: xlatc=76.  !Central latutude of mother domain 
      real, parameter :: xlonc=45.  !central longitude of mother domain

      real :: pi,rpd,rearth,ciy,cjx,ihm,cotrue1,c1,yc,dlon,ypoint,xpoint

      real, intent(in)  :: lon,lat
      real, intent(out) :: x,y
      character(len=*), intent(in) :: grid_id

     pi    =atan2(1.,1.)*4.
     rpd   =pi/180. ! radians per degree
     rearth=6370.949 ! radius of planet, in km
     ciy   =0.5*(1+ydim)
     cjx   =0.5*(1+xdim)

     ihm=1; !Norhtern  hemisphere: ihm=1 else ihm=-1	    
     cotrue1=(ihm*90.)-true1;
     c1=1+cos(rpd*cotrue1)
     yc=-rearth*sin(pi/2. -rpd*xlatc)* c1/(1.+cos( pi/2-rpd*xlatc ));

     dlon=lon-xlonc;




      if (dlon <-180.) dlon=dlon+360.
      if (dlon > 180.) dlon=dlon-360.
      ypoint=-rearth*sin(pi/2.-rpd*lat)* c1/(1.+cos(.5*ihm*pi-rpd*lat))*cos(rpd*dlon)
      xpoint=ihm*rearth*sin(.5*ihm*pi-rpd* lat)*c1/(1.+cos(pi/2-rpd*lat))*sin(rpd*dlon)
      y=(ypoint-yc)/dskmc+ciy
      x=xpoint/dskmc+cjx;
   end subroutine

   subroutine storm_xy2ll(x,y,grid_id,lon,lat)
      implicit none
      real, parameter :: dskmc=12.;      !Gitteravstand i km
      real, parameter :: ydim=230.;      !y dimensions of mother domain
      real, parameter :: xdim=260.;      !x dimensions of mother domain
      real, parameter :: true1=60.;      ! true latitude of polar stereographic grid --- standardized for mm5
      real, parameter :: true2=60.;
      real, parameter :: xlatc=76.0000;  !Central latutude of mother domain 
      real, parameter :: xlonc=45.0000;  ! central longitude of mother domain

      real, intent(out) :: lon,lat
      real, intent(in) :: x,y
      character(len=*), intent(in) :: grid_id

      real :: pi,rpd,rearth,ciy,cjx,ihm,cotrue1,c1,yc,x2,y2

     pi=atan2(1.,1.)*4.; 
     rpd=pi/180.;     ! radians per degree
     rearth=6370.949;! radius of planet, in km
     ciy=0.5*(1+ydim)
     cjx=0.5*(1+xdim) 

     ihm=1; !Norhtern  hemisphere: ihm=1 else ihm=-1	    
     cotrue1=(ihm*90.)-true1;
     c1=1+cos(rpd*cotrue1); 
     yc=-rearth*sin(pi/2 -rpd*xlatc)* c1/(1+cos( pi/2  -  rpd*xlatc )); 

     y2=(y-ciy)*dskmc+yc;
     x2=(x-cjx)*dskmc;
     lat=(pi/2- 2*atan(sqrt(x2**2+y2**2)/(rearth*c1)))/rpd;
     if (abs(x2)<1e-5 .and. abs(y2)<1e-5) then
         lon=xlonc
     else
         lon=xlonc+(atan2(x2,-y2))/rpd
      endif

   end subroutine 

end module mod_storm
