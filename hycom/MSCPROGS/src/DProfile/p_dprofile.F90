program dprofile
   use mod_xc
   use mod_za
   use mod_grid
   use mod_hycomfile_io
   use m_fields_to_plot
   use netcdf
   use m_handle_err
   implicit none
#if defined (IARGC)
   integer*4, external :: iargc
#endif

   character(len=80) :: fnamein, tmparg, ftype
   integer i,j,k,mini,minj, kdm
   real :: rtime,lon,lat,dist, mindist
   real, allocatable :: tmpfld(:,:), pres(:,:,:), tmpfld2(:,:)

   ! Netcdf variable, dimension and file ids
   integer :: ncid, time_dim, kdm_dim, varid_pres, &
      varid(1000)
   character(len=8) :: uname, vname

   type(fields) :: fld(1000)
   integer :: nfld, ifld
   type(hycomfile) :: hfile
   real, external :: spherdist



   ! Two first arguments are position (lon - lat)
   if (iargc()>2) then
      call getarg(1,tmparg) ; read(tmparg,*) lon
      call getarg(2,tmparg) ; read(tmparg,*) lat
   else
      write(6,*) 'Program extracts a profile at a given position.'
      write(6,*) 'You must specify a longitude, latitude pair and'
      write(6,*) 'at least one filename'
      print *
      write(6,*) 'Usage: dprofile [lon] [lat] files...'
      write(6,*) ''
      write(6,*) 'Fields to extract are set in extract - files'
      write(6,*) '( extract.daily, extract.weekly etc etc)'
      stop '(DProfile)'
   end if

   ! Initialize io
   call xcspmd()
   call zaiost()
   call get_grid()

   ! What file is this? (daily? weekly? restart? pak?)
   call getarg(3,fnamein) 
   ftype=getfiletype(trim(fnamein))

   ! Inits file type
   call initHF(hfile,fnamein,trim(ftype))
   kdm=vDim(hfile)
   allocate(tmpfld (idm,jdm))
   allocate(tmpfld2(idm,jdm))
   allocate(pres(idm,jdm,kdm+1))

   ! Get fields to plot
   call fields_to_plot(fld,nfld,hfile,kdm)

   ! Locate closest point on grid
   mindist=1e10
   mini=-1
   minj=-1
   do j=1,jdm
   do i=1,idm
      dist=spherdist(lon,lat,plon(i,j),plat(i,j))
      if (dist< mindist) then
         mini=i
         minj=j
         mindist=dist
      end if
   end do
   end do
   print '(a,2f14.2)','Requested position : ',lon,lat
   print '(a,2f14.2)','got       position : ',plon(mini,minj),plat(mini,minj)

   if (NF90_CREATE('profile.nc',NF90_CLOBBER,ncid) /= NF90_NOERR) then
      print *,'An error occured when opening the netcdf file'
      stop '(norsexclim)'
   end if
   call handle_err(NF90_DEF_DIM(ncid,'time',NF90_UNLIMITED,time_dim))
   call handle_err(NF90_DEF_DIM(ncid,'depth',kdm,kdm_dim))
   call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'longitude',lon))
   call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'latitude',lat))
   call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'actual_longitude',plon (mini,minj)))
   call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'actual_latitude' ,plat (mini,minj)))
   call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'model_depth'     ,depths(mini,minj)))

   ! Define all sought-for scalar variables
   call handle_err(NF90_DEF_VAR(ncid,'pres', NF90_Float,(/kdm_dim,time_dim/),varid_pres))
   do ifld=1,nfld
   if (fld(ifld)%option .and. (.not. fld(ifld)%vecflag)) then
      if (is3DVar(hfile,fld(ifld)%fextract,1)) then
         ! Define variable on first pass
         call handle_err(NF90_DEF_VAR(ncid,fld(ifld)%fextract, &
            NF90_Float,(/kdm_dim,time_dim/),varid(ifld)))
      else
         ! Define variable on first pass
         call handle_err(NF90_DEF_VAR(ncid,fld(ifld)%fextract, &
            NF90_Float,(/time_dim/),varid(ifld)))
      end if
   end if
   end do

   ! Define all sought-for vector variables
   do ifld=1,nfld-1
   if (fld(ifld)%option .and. fld(ifld)%vecflag) then
      uname=fld(ifld  )%fextract
      vname=fld(ifld+1)%fextract
      if (is3DVar(hfile,fld(ifld)%fextract,1)) then
         call handle_err(NF90_DEF_VAR(ncid,trim(uname)//'_east', &
            NF90_Float,(/kdm_dim,time_dim/),varid(ifld  )))
         call handle_err(NF90_DEF_VAR(ncid,trim(vname)//'_north', &
            NF90_Float,(/kdm_dim,time_dim/),varid(ifld+1)))
      else
         call handle_err(NF90_DEF_VAR(ncid,trim(uname)//'_east', &
            NF90_Float,(/time_dim/),varid(ifld )))
         call handle_err(NF90_DEF_VAR(ncid,trim(vname)//'_north', &
            NF90_Float,(/time_dim/),varid(ifld+1)))
      end if
   end if
   end do



   call handle_err(NF90_enddef(ncid))


   ! Cycle files
   do i = 3,iargc()
      call getarg(i,fnamein) 
      call initHF(hfile,fnamein,trim(ftype))
      print '(a)','processing '//trim(fnamein)

      ! Get pressures
      pres(:,:,1)=0.
      do k=1,kdm
         !call HFReadDPField(hfile,tmpfld,idm,jdm,k,1)
         call HFReadDPField_m(hfile,tmpfld,idm,jdm,k,1) ! Gives dp in meters
         where(tmpfld>0.5*huge) tmpfld=0.
         pres(:,:,k+1)=pres(:,:,k)+tmpfld
      end do
      !pres=pres/onem

      ! Put pressure variables
      call handle_err(NF90_PUT_VAR(ncid,varid_pres, &
         pres(mini,minj,2:kdm+1),start=(/1,i-2/)))


      ! Put scalar variables
      do ifld=1,nfld
      if (fld(ifld)%option .and. (.not. fld(ifld)%vecflag)) then
         print '(a)','--Extracting '//fld(ifld)%fextract
         if (is3DVar(hfile,fld(ifld)%fextract,1)) then
            do k=1,kdm
               call HFReadField(hfile,tmpfld,idm,jdm,fld(ifld)%fextract,k,1)
               call handle_err(NF90_PUT_VAR(ncid,varid(ifld), &
                               tmpfld(mini,minj),start=(/k,i-2/)))
            end do
         else
            call HFReadField(hfile,tmpfld,idm,jdm,fld(ifld)%fextract,0,1)
            call handle_err(NF90_PUT_VAR(ncid,varid(ifld), &
                            tmpfld(mini,minj),start=(/i-2/)))
         end if
      end if ! process flag
      end do ! Field loop


      ! Put vector variables
      do ifld=1,nfld-1
      if (fld(ifld)%option .and. fld(ifld)%vecflag) then
         uname=fld(ifld  )%fextract
         vname=fld(ifld+1)%fextract
         print '(a)','--Extracting '//uname//vname
         if (is3DVar(hfile,fld(ifld)%fextract,1)) then
            do k=1,kdm
               call HFReadField(hfile,tmpfld ,idm,jdm,uname,k,1)
               call HFReadField(hfile,tmpfld2,idm,jdm,vname,k,1)
               call rotate(tmpfld,tmpfld2,plat,plon,idm,jdm,'m2l')

               call handle_err(NF90_PUT_VAR(ncid,varid(ifld  ), &
                               tmpfld (mini,minj),start=(/k,i-2/)))
               call handle_err(NF90_PUT_VAR(ncid,varid(ifld+1), &
                               tmpfld2(mini,minj),start=(/k,i-2/)))
            end do
         else
            call HFReadField(hfile,tmpfld ,idm,jdm,uname,0,1)
            call HFReadField(hfile,tmpfld2,idm,jdm,vname,0,1)
            call rotate(tmpfld,tmpfld2,plat,plon,idm,jdm,'m2l')
            call handle_err(NF90_PUT_VAR(ncid,varid(ifld), &
                            tmpfld (mini,minj),start=(/i-2/)))
            call handle_err(NF90_PUT_VAR(ncid,varid(ifld+1), &
                            tmpfld2(mini,minj),start=(/i-2/)))
         end if
      end if ! process flag
      end do ! Field loop

   end do ! File loop i
   call handle_err(NF90_close(ncid))
   print *,'profiles in profile.nc'
end program dprofile
