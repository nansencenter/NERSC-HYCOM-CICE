! Program to convert from various hycom - type files
! to netcdf files with data interploated along a section
! KAL: Added options to read other than pak files

program section_plot
   use mod_xc
   use mod_za
   use mod_grid
   use mod_year_info
   use mod_hycomfile_io
   use m_fields_to_plot
   use m_handle_err
   use netcdf
   implicit none

   character(len=24) :: char24,plot_time,tmparg
   character(len= 8) :: varin,tmpc
   character(len= 3) :: rungen


   integer :: utj,vtj

   integer i,j,ia,ib,ja,jb,n2d,n3d, kmin, kstart
   integer,allocatable :: seccnt(:,:)
   integer ii,jj,kk,ii1,jj1,k,n,nrrec,l
   integer lay1,lay2,ndim,nstep
   integer, parameter :: iversion=1
   integer ivers,iday,iyear,ihour,imonth,iweek
   character(len=5) dpfield
   integer, parameter :: maxch=70
   character(len=5), dimension(maxch) :: ch2d,ch3d
   real :: rtime
   integer :: velrot, klevel
   logical :: test_normal_tangent

   integer :: hycfile,nrmem,kdm
   character(len=80) :: fnamein, fbase,newname, ftype

   type(year_info) :: rti, rtd
   real, allocatable :: twod(:,:), fldzonal(:), areazonal(:), &
      dp(:,:), field(:,:), latzonal(:), field2(:,:), fldzonal2(:), &
      pres(:,:,:), field3d(:,:,:), field3d2(:,:,:), fldsum(:,:), dpsum(:,:), &
      fldsum2(:,:)
   logical, allocatable:: zonalmask(:,:), listmask(:,:)
   real :: minlat, maxlat

   real :: urot,vrot,dlon,dlat,utmp,vtmp,theta_up,theta_vp,theta_sp
   real :: radinv,crossprod(3)
   integer :: i0,j0,ipnt,isec,iarg
   logical :: appendfile, isvec

   integer :: ilat,npts,nlat
   real :: lwlat,uplat,cntlat, mindx, cnt
   type(hycomfile) :: hfile
   type(fields) :: fld(1000)
   integer :: nfld, ifld

   real, external :: spherdist
   real, parameter :: undefr4=real(undef,kind=4)
   integer :: ncid, latdim_id, vdim_id, rdim_id, dimms1D(1),  &
      dimms2D(2), dimms3D(3), varid,var(1000)

   real, dimension(:), allocatable :: deep(:)
   integer :: ndeep, ideep
   real, parameter  ::dstep=50.
   real, parameter :: dxfac=3


   ndeep= 5000./dstep
   allocate(deep(ndeep))
   do k=1,ndeep
      deep(k)=k*dstep
   end do




   radinv=pi/180.


   call getarg(1,fnamein)

   ! What file is this? (daily? weekly? restart? pak?)
   ! TODO: force filetype option goes here
   ftype=getfiletype(trim(fnamein))

   ! Inits file type
   call initHF(hfile,trim(fnamein),trim(ftype))
   kdm=vdim(hfile)
   rtime=hfile%fyear

   ! Initialize IO for .ab files
   CALL XCSPMD()  
   CALL ZAIOST()
   call get_grid()


   ! Retrieve fields to extract
   call fields_to_plot(fld,nfld,hfile,kdm)

   ! Read layer interfaces
   allocate(dp   (idm,jdm))
   allocate(dpsum (idm,jdm))
   allocate(fldsum(idm,jdm))
   allocate(fldsum2(idm,jdm))
   allocate(zonalmask(idm,jdm))
   allocate(listmask (idm,jdm))
   allocate(dp   (idm,jdm))
   allocate(pres (idm,jdm,kdm+1))
   allocate(field (idm,jdm))
   allocate(field2(idm,jdm))
   allocate(field3d (idm,jdm,kdm))
   allocate(field3d2(idm,jdm,kdm))

   ! min grid size
   mindx=min(minval(scuy),minval(scvx))*dxfac ! factor subject to tuning

   ! Convert min dx to latitude increments
   dlat=mindx/spherdist(0.,0.,0.,1.)
   minlat=minval(plat)
   maxlat=maxval(plat)
   ! Security check 1
   if (mindx<1000.) then
      print *,'Too low min grid spacing :',mindx
      stop
   ! KAL - Added this test because some older regional files can have 1e30 for max lat
   else if (abs(minlat) > 90. .or. abs(maxlat) > 90) then
      print *,'Security check, |lat| > 90  ... Abandon ship'
      stop
   end if

   ! Allocate "zonal storage"
   nlat=(maxlat-minlat)/dlat
   allocate(fldzonal (nlat))
   allocate(fldzonal2(nlat))
   allocate(latzonal(nlat+1))
   do l=1,nlat+1
      latzonal(l)=minlat + (l-1)*dlat
   end do


   ! Set up netcdf file
   if (NF90_create('tmpzonal.nc',NF90_CLOBBER,ncid) /= NF90_NOERR) then
      print *,'An error occured when opening the netcdf file'
      stop '(restart2netcdf)'
   end if



   ! Define dimensions
   call handle_err(NF90_DEF_DIM(ncid,'latitude',nlat          ,latdim_id))
   call handle_err(NF90_DEF_DIM(ncid,'depth',ndeep        ,vdim_id))
   call handle_err(NF90_DEF_DIM(ncid,'rdim',nf90_unlimited,rdim_id))
   dimms1D=(/latdim_id/)
   dimms2D=(/latdim_id,rdim_id/)
   dimms3D=(/latdim_id,vdim_id,rdim_id/)

   ! Define some basic vars
   do l=1,nlat
      fldzonal(l) = sum( depths, mask=plat>latzonal(l) .and.  plat<latzonal(l+1) .and. depths>.1)
      cnt=count(mask=plat>latzonal(l) .and.  plat<latzonal(l+1) .and. depths>.1)
      fldzonal(l)= fldzonal(l)/cnt
   end do
   call handle_err(NF90_DEF_VAR(ncid,'model_depth',NF90_Float,latdim_id,varid))
   call handle_err(NF90_PUT_ATT(ncid,varid,'_FillValue',real(undefr4,kind=4)))
   call handle_err(NF90_PUT_ATT(ncid,varid,'missing_value',undefr4))
   call handle_err(NF90_ENDDEF(ncid))
   call handle_err(NF90_PUT_VAR(ncid,varid,fldzonal))

   call handle_err(NF90_REDEF(ncid))
   call handle_err(NF90_DEF_VAR(ncid,'latitude',NF90_Float,latdim_id,varid))
   call handle_err(NF90_ENDDEF(ncid))
   call handle_err(NF90_PUT_VAR(ncid,varid,0.5*(latzonal(1:nlat)+latzonal(2:nlat+1))))

   call handle_err(NF90_REDEF(ncid))
   call handle_err(NF90_DEF_VAR(ncid,'depth',NF90_Float,vdim_id,varid))
   call handle_err(NF90_ENDDEF(ncid))
   call handle_err(NF90_PUT_VAR(ncid,varid,deep))

   !call handle_err(NF90_CLOSE(ncid))
   !stop


   pres=0.
   do k=1,kdm
      !call HFReadDPField(hfile,dp,idm,jdm,k,1)
      call HFReadDPField_m(hfile,dp,idm,jdm,k,1) ! Always returns dp in meters
      where(dp>0.5*huge) dp=0.
      pres(:,:,k+1)=pres(:,:,k)+dp
   end do
   pres=pres/onem

   ! Define all generic variables
   call handle_err(NF90_REDEF(ncid))
   do ifld=1,nfld
      if (is3DVar(hfile,fld(ifld)%fextract,1)) then
         call handle_err(NF90_DEF_VAR(ncid,fld(ifld)%fextract,NF90_Float,dimms3D,var(ifld)))
         call handle_err(NF90_PUT_ATT(ncid,var(ifld),'_FillValue',real(undef,kind=4)))
      else
         call handle_err(NF90_DEF_VAR(ncid,fld(ifld)%fextract,NF90_Float,dimms2D,var(ifld)))
         call handle_err(NF90_PUT_ATT(ncid,var(ifld),'_FillValue',real(undef,kind=4)))
      end if
   end do
   call handle_err(NF90_ENDDEF(ncid))



   ! Go through fields to extract
   !TODO: Rotation doesnt change variable names
   !TODO: Weight with grid size
   do ifld=1,nfld
      if (is3DVar(hfile,fld(ifld)%fextract,1)) then

         ! Special case (and rules) for vectors -
         ! 1) No scalar variables should begin with u,v,taux or tauy
         ! 2) v-component must come immediately after u-component
         ! See also checks in  m_fields_to_plot
         isvec = fld(ifld)%fextract(1:1)=='u' .or. fld(ifld)%fextract(1:4)=='taux'
         if (fld(ifld)%fextract(1:1)/='v' .and. fld(ifld)%fextract(1:4)/='tauy') then
            write(6,'(a)') '3D Var: '//trim(fld(ifld)%fextract)

            ! Read 3D field
            call HFReadField3D(hfile,field3d ,idm,jdm,kdm,fld(ifld  )%fextract,1)
            if (isvec)  then
               call HFReadField3D(hfile,field3d2 ,idm,jdm,kdm,fld(ifld+1)%fextract,1)
               do k=1,kdm
                  call rotate(field3d(:,:,k),field3d2(:,:,k),plat,plon,idm,jdm,'m2l')
               end do
            end if


            ! for each deep, integrate between prev interface and next
            ! TODO: sloooow
            kstart=1
            do ideep=1,ndeep
               !print *,ideep
               dpsum=0.
               fldsum=0.
               fldsum2=0.
               do k=kstart,kdm
                  if (ideep==1) then
                     dp=max(0.,min(pres(:,:,k+1),deep(ideep))-pres(:,:,k))
                  else
                     dp=max(0.,min(pres(:,:,k+1),deep(ideep))-max(pres(:,:,k),deep(ideep-1)))
                  end if

                  ! Flag 1
                  listmask=pres(:,:,k+1)< deep(ideep) .or. pres(:,:,kdm+1) < deep(ideep)
                  if (count(listmask)==idm*jdm) kmin=k

                  !! Flag 2
                  !listmask=pres(:,:,k)> deep(ideep) .or. pres(:,:,kdm+1) < deep(ideep)
                  !if (count(listmask)==idm*jdm) exit

                  dpsum=dpsum+dp
                  fldsum=fldsum+dp*field3d(:,:,k)
                  if (isvec) fldsum2=fldsum2+dp*field3d2(:,:,k)
               end do
               print *,kmin
               kstart=kmin

               where (dpsum>10.) 
                  fldsum=fldsum/dpsum
               elsewhere
                  fldsum=undef
               endwhere

               if (isvec) then
                  where (dpsum>10.) 
                     fldsum2=fldsum2/dpsum
                  elsewhere
                     fldsum2=undef
                  endwhere
               end if
                     

               ! Extract in meridional bands. Straight-forward, and inefficient way ! - fix if it bugs you !
               do l=1,nlat
                  zonalmask=plat>latzonal(l) .and.  plat<latzonal(l+1) .and. depths>.1 .and. fldsum/=undef .and.fldsum2/=undef
                  fldzonal (l) = sum( fldsum , mask=zonalmask)
                  if (isvec) fldzonal2(l) = sum( fldsum2, mask=zonalmask)
                  cnt=count(mask=zonalmask)
                  if (cnt>0) then
                     fldzonal (l)= fldzonal(l)/cnt
                     if (isvec) fldzonal2(l)= fldzonal2(l)/cnt
                  else
                     fldzonal (l)= undef
                     if (isvec) fldzonal2(l)= undef
                  end if
                  !end if
               end do
               call handle_err(NF90_PUT_VAR(ncid,var(ifld  ),fldzonal ,start=(/1,ideep,1/)))
               if (isvec) call handle_err(NF90_PUT_VAR(ncid,var(ifld+1),fldzonal2 ,start=(/1,ideep,1/)))
            end do


!            do k=1,kdm
!               call HFReadField(hfile,field ,idm,jdm,fld(ifld  )%fextract,k,1)
!               where (field >0.5*huge) field=0.
!               if (isvec) then
!                  call HFReadField(hfile,field2,idm,jdm,fld(ifld+1)%fextract,k,1)
!                  where (field2>0.5*huge) field2=0.
!                  call rotate(field,field2,plat,plon,idm,jdm,'m2l')
!                  !write(6,'(a)') '3D Var: '//trim(fld(ifld+1)%fextract)//' rotate'
!               end if
!
!               ! Extract in meridional bands. Straight-forward, and inefficient way ! - fix if it bugs you !
!               do l=1,nlat
!                  fldzonal (l) = sum( field , mask=plat>latzonal(l) .and.  plat<latzonal(l+1) .and. depths>.1)
!                  cnt=count(mask=plat>latzonal(l) .and.  plat<latzonal(l+1) .and. depths>.1)
!                  fldzonal (l)= fldzonal(l)/cnt
!                  if (isvec) then
!                     fldzonal2(l) = sum( field2, mask=plat>latzonal(l) .and.  plat<latzonal(l+1) .and. depths>.1)
!                     fldzonal2(l)= fldzonal2(l)/cnt
!                  end if
!               end do
!               call handle_err(NF90_PUT_VAR(ncid,var(ifld  ),fldzonal ,start=(/1,k,1/)))
!               if(isveC) call handle_err(NF90_PUT_VAR(ncid,var(ifld+1),fldzonal2,start=(/1,k,1/)))
!            end do
         end if

      ! 2D case
      else

         ! Special case (and rules) for vectors -
         ! 1) No scalar variables should begin with u,v,taux or tauy
         ! 2) v-component must come immediately after u-component
         ! See also checks in  m_fields_to_plot
         isvec = fld(ifld)%fextract(1:1)=='u' .or. fld(ifld)%fextract(1:4)=='taux'
         if (fld(ifld)%fextract(1:1)/='v' .and. fld(ifld)%fextract(1:4)/='tauy') then
            write(6,'(a)') '2D Var: '//trim(fld(ifld)%fextract)
            call HFReadField(hfile,field ,idm,jdm,fld(ifld  )%fextract,0,1)
            where (field>0.5*huge) field=0.
            if (isvec) then
               call HFReadField(hfile,field2,idm,jdm,fld(ifld+1)%fextract,0,1)
               where (field2>0.5*huge) field2=0.
               call rotate(field,field2,plat,plon,idm,jdm,'m2l')
            end if

            ! Extract in meridional bands. Straight-forward, and inefficient way
            ! - fix if it bugs you !
            do l=1,nlat
               fldzonal (l) = sum( field , mask=plat>latzonal(l) .and.  plat<latzonal(l+1) .and. depths>.1)
               cnt=count(mask=plat>latzonal(l) .and.  plat<latzonal(l+1) .and. depths>.1)
               fldzonal (l)= fldzonal (l)/cnt
               if (isvec) then
                  fldzonal2(l) = sum( field2, mask=plat>latzonal(l) .and.  plat<latzonal(l+1) .and. depths>.1)
                  fldzonal2(l)= fldzonal2(l)/cnt
                  !write(6,'(a)') '2D Var: '//trim(fld(ifld+1)%fextract)//' rotate'
               end if
            end do
            call handle_err(NF90_PUT_VAR(ncid,var(ifld  ),fldzonal ,start=(/1,1/)))
            if (isvec) call handle_err(NF90_PUT_VAR(ncid,var(ifld+1),fldzonal2,start=(/1,1/)))

         end if
      end if
   end do
   call handle_err(NF90_CLOSe(ncid))



end program section_plot
