      program check_nest
      use netcdf
      use m_get_nest_record
      use m_get_nest_info
      implicit none

      character(len= *),parameter :: fnestloc='nestloc.uf'
      character(len=80)           :: fndepths,nestingfile,matfil, ncfile
      character(len= 7)           :: tag7
      character(len=20)           :: ctitle

      ! Dimensions of inner grid
      integer :: idm, jdm, inest, &
         iidm, jjdm

      ! Arrays holding data on inner grid (where we need nesting
      ! conditions)
      real*8, allocatable, dimension(:,:) :: &
         inner_ior8 ,mat_iofld
      real, allocatable, dimension(:,:) :: tmpfld
      real, allocatable, dimension(:,:) :: &
         inner_lon, inner_lat, inner_depths
      integer, allocatable, dimension(:,:) :: &
         inner_ipiv, inner_jpiv

      real, parameter  :: onem=9806.
      integer :: ixx,ix,jx

      integer :: irec,klevel, inum_offset,ivar, &
         kdm, intvar, nrec_offset, num_offset, firstrec, lastrec
      logical :: ex, isvelocity,ass

      integer :: i1dim,i2dim,j1dim,j2dim, ibnd, itime
      integer :: i1,i2,j1,j2,indbnd,imax,iday,iyear,i,j,k,l,k2
      character(len=2) :: cbnd
      character(len=4) :: cll
      character(len=3) :: cvar
      integer :: ncid, varid, idmid, jdmid, rdimid, ierr, kdmid
      real, parameter :: undef=-1e14

      inquire(file=fnestloc,exist=ex)
      if (.not.ex) then 
         print *,fnestloc//' does not exist!'
         stop '(nest_offline)'
      end if


      ! Read nesting positions for the internal grid
      open(10,file=fnestloc,form='unformatted',status='old')

      ! Read grid dimensions
      read(10)idm,jdm,inest
      iidm=idm-inest+1
      jjdm=jdm-inest+1

      write(6,'(a,5i5)') &
         fnestloc//'  has dimensions: ', idm,jdm,&
         iidm,jjdm,inest

      ! Allocate grid for inner model
      allocate(inner_lon(idm,jdm))
      allocate(inner_lat(idm,jdm))
      allocate(inner_ior8(idm,jdm))
      allocate(tmpfld    (idm,jdm))

      ! Read grid
      read(10) inner_ior8 ; inner_lon=inner_ior8
      read(10) inner_ior8 ; inner_lat=inner_ior8
      close(10)

      write(6,'(a)')   fnestloc//'  is read'

      ! Allocate temporary arrays for interpolating to inner grid
      ! + depth matrix
      allocate(inner_depths(idm,jdm))


      ! Read the depth matrix of the local grid 
! --- tag7 has 11 chars to accomodate for huge grids in future
      if (idm>999 .or. jdm > 999) then
         write(tag7,'(i5.5,a,i5.5)')idm,'x',jdm
      else
         write(tag7,'(i3.3,a,i3.3)')idm,'x',jdm
      end if

      fndepths='ndepths'//trim(tag7)//'.uf'

      inquire(file=trim(fndepths),  exist=ex)
      if (.not.ex) then
          write(6,'(a)') &
         'nesting depths file for local grid does not exist:', &
                 'ndepths'//trim(tag7)//'.uf'
         stop '(nest_offline)'
      else
         open(10,file=trim(fndepths),form='unformatted',status='old')
         read(10) inner_ior8
         close(10)
         inner_depths=inner_ior8
      endif
         
      ! Extend local grid to boundary
      where (inner_depths(2,:)>0.)  &
         inner_depths(1,:)=inner_depths(2,:)
      where (inner_depths(idm-1,:)>0.)  &
         inner_depths(idm,:)=inner_depths(idm-1,:)
      where (inner_depths(:,2)>0.)  &
         inner_depths(:,1)=inner_depths(:,2)
      where (inner_depths(:,jdm-1)>0.) &
         inner_depths(:,jdm)=inner_depths(:,jdm-1)


      ! Get date to process
      print *,'Input year and julian day: '
      read (*,*) iyear, iday

      write(nestingfile,'(a,i4.4,a,i3.3)') 'Nest/nest_',iyear,'_',iday

      ! Get kdm and offsets from header file
      call get_nest_info(trim(nestingfile),kdm,nrec_offset,num_offset)
      print *,'From get_nest_info:'
      print *,'kdm(number of layers)          :',kdm
      print *,'nrec_offset(records per time  ):',nrec_offset
      print *,'num_offset (time dumps in file):',num_offset

      write(ncfile,'(a,i4.4,a,i3.3,a)') 'nest_',iyear,'_',iday,'.nc'
      if (NF90_CREATE(trim(ncfile),NF90_CLOBBER,ncid) /= NF90_NOERR) then
         print *,'An error occured when opening the netcdf file'
         stop '(obsstats)'
      end if
      ierr=NF90_DEF_DIM(ncid,'idm',idm,idmid)
      ierr=NF90_DEF_DIM(ncid,'jdm',jdm,jdmid)
      ierr=NF90_DEF_DIM(ncid,'kdm',kdm,kdmid)
      ierr=NF90_DEF_DIM(ncid,'rdm',NF90_UNLIMITED,rdimid)

      ierr=NF90_DEF_VAR(ncid,'longitude',NF90_Float,(/idmid,jdmid/),varid)
      ierr=NF90_ENDDEF(ncid)
      ierr=NF90_PUT_VAR(ncid,varid,inner_lon)
      ierr=NF90_REDEF(ncid)

      ierr=NF90_DEF_VAR(ncid,'latitude',NF90_Float,(/idmid,jdmid/),varid)
      ierr=NF90_ENDDEF(ncid)
      ierr=NF90_PUT_VAR(ncid,varid,inner_lat)
      ierr=NF90_REDEF(ncid)

      ierr=NF90_DEF_VAR(ncid,'depth',NF90_Float,(/idmid,jdmid/),varid)
      ierr=NF90_ENDDEF(ncid)
      ierr=NF90_PUT_VAR(ncid,varid,inner_depths)
      ierr=NF90_REDEF(ncid)


         
      print *,'interface ...'
      ierr=NF90_DEF_VAR(ncid,'interface',NF90_Float,(/idmid,jdmid,kdmid,rdimid/),varid)
      ierr=NF90_PUT_ATT(ncid,varid,'_FillValue',undef)
      ierr=NF90_PUT_ATT(ncid,varid,'missing_value',undef)
      ierr=NF90_ENDDEF(ncid)
      do itime=1,num_offset
         do k=1,kdm
            call read_nest('INT',k,nestingfile,tmpfld,idm,jdm,inest,itime)
            ierr=NF90_PUT_VAR(ncid,varid,tmpfld,start=(/1,1,k,itime/))
         end do
      end do
      ierr=NF90_REDEF(ncid)

      print *,'temperature ...'
      ierr=NF90_DEF_VAR(ncid,'temperature',NF90_Float,(/idmid,jdmid,kdmid,rdimid/),varid)
      ierr=NF90_PUT_ATT(ncid,varid,'_FillValue',undef)
      ierr=NF90_PUT_ATT(ncid,varid,'missing_value',undef)
      ierr=NF90_ENDDEF(ncid)
      do itime=1,num_offset
         do k=1,kdm
            call read_nest('TEM',k,nestingfile,tmpfld,idm,jdm,inest,itime)
            ierr=NF90_PUT_VAR(ncid,varid,tmpfld,start=(/1,1,k,itime/))
         end do
      end do
      ierr=NF90_REDEF(ncid)

      print *,'salinity ...'
      ierr=NF90_DEF_VAR(ncid,'salinity',NF90_Float,(/idmid,jdmid,kdmid,rdimid/),varid)
      ierr=NF90_PUT_ATT(ncid,varid,'_FillValue',undef)
      ierr=NF90_PUT_ATT(ncid,varid,'missing_value',undef)
      ierr=NF90_ENDDEF(ncid)
      do itime=1,num_offset
         do k=1,kdm
            call read_nest('SAL',k,nestingfile,tmpfld,idm,jdm,inest,itime)
            ierr=NF90_PUT_VAR(ncid,varid,tmpfld,start=(/1,1,k,itime/))
         end do
      end do
      ierr=NF90_REDEF(ncid)

      print *,'u baroclinic (eastward) ...'
      ierr=NF90_DEF_VAR(ncid,'ueast',NF90_Float,(/idmid,jdmid,kdmid,rdimid/),varid)
      ierr=NF90_PUT_ATT(ncid,varid,'_FillValue',undef)
      ierr=NF90_PUT_ATT(ncid,varid,'missing_value',undef)
      ierr=NF90_ENDDEF(ncid)
      do itime=1,num_offset
         do k=1,kdm
            call read_nest('UT',k,nestingfile,tmpfld,idm,jdm,inest,itime)
            ierr=NF90_PUT_VAR(ncid,varid,tmpfld,start=(/1,1,k,itime/))
         end do
      end do
      ierr=NF90_REDEF(ncid)

      print *,'v baroclinic (northward) ...'
      ierr=NF90_DEF_VAR(ncid,'vnorth',NF90_Float,(/idmid,jdmid,kdmid,rdimid/),varid)
      ierr=NF90_PUT_ATT(ncid,varid,'_FillValue',undef)
      ierr=NF90_PUT_ATT(ncid,varid,'missing_value',undef)
      ierr=NF90_ENDDEF(ncid)
      do itime=1,num_offset
         do k=1,kdm
            call read_nest('VT',k,nestingfile,tmpfld,idm,jdm,inest,itime)
            ierr=NF90_PUT_VAR(ncid,varid,tmpfld,start=(/1,1,k,itime/))
         end do
      end do
      ierr=NF90_REDEF(ncid)

      print *,'u barotropic (eastward) ...'
      ierr=NF90_DEF_VAR(ncid,'ubeast',NF90_Float,(/idmid,jdmid,rdimid/),varid)
      ierr=NF90_PUT_ATT(ncid,varid,'_FillValue',undef)
      ierr=NF90_PUT_ATT(ncid,varid,'missing_value',undef)
      ierr=NF90_ENDDEF(ncid)
      do itime=1,num_offset
         call read_nest('UB',0,nestingfile,tmpfld,idm,jdm,inest,itime)
         ierr=NF90_PUT_VAR(ncid,varid,tmpfld,start=(/1,1,itime/))
      end do
      ierr=NF90_REDEF(ncid)

      print *,'v barotropic (northward) ...'
      ierr=NF90_DEF_VAR(ncid,'vbnorth',NF90_Float,(/idmid,jdmid,rdimid/),varid)
      ierr=NF90_PUT_ATT(ncid,varid,'_FillValue',undef)
      ierr=NF90_PUT_ATT(ncid,varid,'missing_value',undef)
      ierr=NF90_ENDDEF(ncid)
      do itime=1,num_offset
         call read_nest('VB',0,nestingfile,tmpfld,idm,jdm,inest,itime)
         ierr=NF90_PUT_VAR(ncid,varid,tmpfld,start=(/1,1,itime/))
      end do
      ierr=NF90_REDEF(ncid)


      print *,'fice  ...'
      ierr=NF90_DEF_VAR(ncid,'hice',NF90_Float,(/idmid,jdmid,rdimid/),varid)
      ierr=NF90_PUT_ATT(ncid,varid,'_FillValue',undef)
      ierr=NF90_PUT_ATT(ncid,varid,'missing_value',undef)
      ierr=NF90_ENDDEF(ncid)
      do itime=1,num_offset
         call read_nest('HI',0,nestingfile,tmpfld,idm,jdm,inest,itime)
         ierr=NF90_PUT_VAR(ncid,varid,tmpfld,start=(/1,1,itime/))
      end do
      ierr=NF90_REDEF(ncid)

      print *,'hice  ...'
      ierr=NF90_DEF_VAR(ncid,'hice',NF90_Float,(/idmid,jdmid,rdimid/),varid)
      ierr=NF90_PUT_ATT(ncid,varid,'_FillValue',undef)
      ierr=NF90_PUT_ATT(ncid,varid,'missing_value',undef)
      ierr=NF90_ENDDEF(ncid)
      do itime=1,num_offset
         call read_nest('HI',0,nestingfile,tmpfld,idm,jdm,inest,itime)
         ierr=NF90_PUT_VAR(ncid,varid,tmpfld,start=(/1,1,itime/))
      end do
      ierr=NF90_REDEF(ncid)

!      print *,'hsnw  ...'
!      ierr=NF90_DEF_VAR(ncid,'hsnw',NF90_Float,(/idmid,jdmid,rdimid/),varid)
!      ierr=NF90_PUT_ATT(ncid,varid,'_FillValue',undef)
!      ierr=NF90_PUT_ATT(ncid,varid,'missing_value',undef)
!      ierr=NF90_ENDDEF(ncid)
!      do itime=1,num_offset
!         call read_nest('HS',0,nestingfile,tmpfld,idm,jdm,inest,itime)
!         ierr=NF90_PUT_VAR(ncid,varid,tmpfld,start=(/1,1,itime/))
!      end do
!      ierr=NF90_REDEF(ncid)

      ierr=NF90_CLOSE(ncid)
         






      stop '(normal)'

     contains

        subroutine read_nest(cvar,k,nestingbase,fld,nx,ny,inest,itime)
        implicit none
        integer         , intent(in) :: inest, nx, ny, itime
        character(len=*), intent(in) :: cvar,nestingbase
        integer         , intent(in) :: k
        real, intent(out)            :: fld(nx,ny)

        integer :: i1,i2,j1,j2, irec
        real*4 :: inner_ior4(nx, ny)


         firstrec=(itime-1)*nrec_offset+1
         lastrec =(itime  )*nrec_offset
         fld=undef

         !print *,cvar, k

         call get_nest_record(cvar,k,trim(nestingbase),irec,firstrec,lastrec)

         if (irec<0) then 
            print '(a,i3,a)','Variable '//cvar,k,' not found'
            fld=undef
            return
         end if


         do ibnd=1,4
            if     (ibnd==1)  then
               cbnd='i1'
               i1=1    ; i2=inest
               j1=1    ; j2=ny
            elseif (ibnd==2)  then
               cbnd='j1'
               i1=1    ; i2=nx
               j1=1    ; j2=inest
            elseif (ibnd==3)  then
               cbnd='ii'
               i1=nx-inest+1   ; i2=nx
               j1=1    ; j2=ny
            elseif (ibnd==4)  then
               cbnd='jj'
               i1=1    ; i2=nx
               j1=ny-inest+1   ; j2=ny
            end if

            ! Open nesting file
!            write(6,'(a,a,a)') 'Reading records ',  &
!               ' from file=',trim(nestingfile)//'_'//cbnd
            inquire(iolength=j)inner_ior4(i1:i2,j1:j2)
            open(10,file=trim(nestingfile)//'_'//cbnd,form='unformatted', &

              access='direct',recl=j,status='old')
            read(10,rec=irec) inner_ior4(i1:i2,j1:j2)
            fld(i1:i2,j1:j2)=inner_ior4(i1:i2,j1:j2)
            close(10)
         end do
         where(inner_depths<.1) fld=undef
         end subroutine


!
!      ! Tecplot IO
!      imax=max(i2-i1,j2-j1)+1
!      open(66,file=trim(nestingfile)//'_'//cbnd//'_test.tec',status='replace')
!      write(66,'(''TITLE= "Nested fields test"'')')
!      write(66,'(a)')'VARIABLES="i" "k" "lon" "lat" "depth" "saln" "temp"'// &
!                                '"u(baclin)" "v(baclin)" "u(total)" "v(total)"'
!      do l=1,inest
!         if (cbnd(1:1)=='i') then
!            i1=i1dim+(l-1)
!            i2=i1dim+(l-1)
!            j1=1
!            j2=jdm
!            imax=jdm
!            write(ctitle,'(a,i4.4)') 'nestbnd',i1
!            write(cll,'(i4.4)') i1
!         elseif (cbnd(1:1)=='j') then
!            j1=j1dim+(l-1)
!            j2=j1dim+(l-1)
!            i1=1
!            i2=idm
!            imax=idm
!            write(ctitle,'(a,i4.4)') 'nestbnd',j1
!            write(cll,'(i4.4)') j1
!         end if
!
!         if (l==1) then
!            WRITE(66,'(a,i6,a,i6)') 'ZONE T='//trim(ctitle)//' F=BLOCK  I=',imax,' J=',kdm
!            WRITE(66,'(30i6)')((i,i=1,imax),j=1,kdm)
!            WRITE(66,'(30i6)')((j,i=1,imax),j=1,kdm)
!            WRITE(66,'(10e15.6)')(((inner_lon(i,j),i=i1,i2),j=j1,j2),k2=1,kdm)
!            WRITE(66,'(10e15.6)')(((inner_lat(i,j),i=i1,i2),j=j1,j2),k2=1,kdm)
!         else if (l>1)  then
!            WRITE(66,'(a,i6,a,i6)') 'ZONE T='//trim(ctitle)//' F=BLOCK  I=',imax,' J=',kdm
!            WRITE(66,'(a)') 'D=(1,2,3,4)'
!         end if
!         print *,'processing '//trim(ctitle)
!
!
!
!
!#if defined(MATLAB)
!         if (.not. allocated (mat_iofld)) allocate(mat_iofld(imax,kdm))
!            ! Get lon -- matlab only
!            do k=1,kdm
!               do jx=j1,j2
!               do ix=i1,i2
!                  ixx=max(ix-i1,jx-j1)+1
!                  mat_iofld(ixx,k)=inner_lon(ix,jx)
!               end do
!               end do
!            end do
!            !print *,irec,maxval(inner_ior4(i1:i2,j1:j2)),i1dim,i2dim,j1dim,j2dim
!            pa1=mxCreateDoubleMatrix(int(imax,kind=IKIND),int(kdm,kind=IKIND),int(0,kind=IKIND))
!            call mxCopyReal8ToPtr(mat_iofld,mxGetPr(pa1),int(imax*kdm,kind=IKIND))
!            iret=matPutVariable(mp, 'longitude'//cll, pa1)
!
!            ! Get lat -- matlab only
!            do k=1,kdm
!               do jx=j1,j2
!               do ix=i1,i2
!                  ixx=max(ix-i1,jx-j1)+1
!                  mat_iofld(ixx,k)=inner_lat(ix,jx)
!               end do
!               end do
!            end do
!            !print *,irec,maxval(inner_ior4(i1:i2,j1:j2)),i1dim,i2dim,j1dim,j2dim
!            !pa1=mxCreateNumericMatrix(imax,kdm,mxClassIDFromClassName('double'),0)
!            call mxCopyReal8ToPtr(mat_iofld,mxGetPr(pa1),int(imax*kdm,kind=IKIND))
!            iret=matPutVariable(mp, 'latitude'//cll, pa1)
!
!         if (l==1) then
!            ! Get lat -- grid distance 
!            do k=1,kdm
!               do jx=j1,j2
!               do ix=i1,i2
!                  ixx=max(ix-i1,jx-j1)+1
!                  mat_iofld(ixx,k)=ixx
!               end do
!               end do
!            end do
!            !print *,irec,maxval(inner_ior4(i1:i2,j1:j2)),i1dim,i2dim,j1dim,j2dim
!            !pa1=mxCreateNumericMatrix(imax,kdm,mxClassIDFromClassName('double'),0)
!            call mxCopyReal8ToPtr(mat_iofld,mxGetPr(pa1),int(imax*kdm,kind=IKIND))
!            iret=matPutVariable(mp, 'grid_distance', pa1)
!         end if
!#endif /*MATLAB*/
!
!
!
!
!         ! Get depths
!         do k=1,kdm
!            call get_nest_record('INT',k,nestingfile,irec,firstrec,lastrec)
!            read(10,rec=irec) inner_ior4(i1dim:i2dim,j1dim:j2dim)
!            where(inner_depths<.1) inner_ior4=0.
!            WRITE(66,'(10e15.6)')((inner_ior4(i,j)/onem,i=i1,i2),j=j1,j2)
!
!            do jx=j1,j2
!            do ix=i1,i2
!               ixx=max(ix-i1,jx-j1)+1
!               mat_iofld(ixx,k)=inner_ior4(ix,jx)
!            end do
!            end do
!         end do
!#if defined(MATLAB)
!         !print *,irec,maxval(inner_ior4(i1:i2,j1:j2)),i1dim,i2dim,j1dim,j2dim
!         call mxCopyReal8ToPtr(mat_iofld,mxGetPr(pa1),int(imax*kdm,kind=IKIND))
!         iret=matPutVariable(mp, 'interface'//cll, pa1)
!#endif /*MATLAB*/
!
!
!         ! Get variable
!         do k=1,kdm
!            call get_nest_record('SAL',k,nestingfile,irec,firstrec,lastrec)
!            read(10,rec=irec) inner_ior4(i1dim:i2dim,j1dim:j2dim)
!            WRITE(66,'(10e15.6)')((inner_ior4(i,j),i=i1,i2),j=j1,j2)
!            !print *,irec,maxval(inner_ior4(i1dim:i2dim,j1dim:j2dim)),i1dim,i2dim,j1dim,j2dim
!
!            do jx=j1,j2
!            do ix=i1,i2
!               ixx=max(ix-i1,jx-j1)+1
!               mat_iofld(ixx,k)=inner_ior4(ix,jx)
!            end do
!            end do
!         end do
!#if defined(MATLAB)
!         call mxCopyReal8ToPtr(mat_iofld,mxGetPr(pa1),int(imax*kdm,kind=IKIND))
!         iret=matPutVariable(mp, 'saln'//cll, pa1)
!#endif /*MATLAB*/
!
!         ! Get variable
!         do k=1,kdm
!            call get_nest_record('TEM',k,nestingfile,irec,firstrec,lastrec)
!            read(10,rec=irec) inner_ior4(i1dim:i2dim,j1dim:j2dim)
!            WRITE(66,'(10e15.6)')((inner_ior4(i,j),i=i1,i2),j=j1,j2)
!            !print *,irec,maxval(inner_ior4(i1dim:i2dim,j1dim:j2dim)),i1dim,i2dim,j1dim,j2dim
!            do jx=j1,j2
!            do ix=i1,i2
!               ixx=max(ix-i1,jx-j1)+1
!               mat_iofld(ixx,k)=inner_ior4(ix,jx)
!            end do
!            end do
!         end do
!#if defined(MATLAB)
!         call mxCopyReal8ToPtr(mat_iofld,mxGetPr(pa1),int(imax*kdm,IKIND))
!         iret=matPutVariable(mp, 'temp'//cll, pa1)
!#endif /*MATLAB*/
!
!         ! Get variable
!         do k=1,kdm
!            call get_nest_record('UT ',k,nestingfile,irec,firstrec,lastrec)
!            read(10,rec=irec) inner_ior4(i1dim:i2dim,j1dim:j2dim)
!            WRITE(66,'(10e15.6)')((inner_ior4(i,j),i=i1,i2),j=j1,j2)
!            !!print *,irec,maxval(inner_ior4(i1dim:i2dim,j1dim:j2dim)),i1dim,i2dim,j1dim,j2dim
!            do jx=j1,j2
!            do ix=i1,i2
!               ixx=max(ix-i1,jx-j1)+1
!               mat_iofld(ixx,k)=inner_ior4(ix,jx)
!            end do
!            end do
!         end do
!#if defined(MATLAB)
!         call mxCopyReal8ToPtr(mat_iofld,mxGetPr(pa1),int(imax*kdm,kind=IKIND))
!         iret=matPutVariable(mp, 'ubaclin'//cll, pa1)
!#endif /*MATLAB*/
!
!         ! Get variable
!         do k=1,kdm
!            call get_nest_record('VT ',k,nestingfile,irec,firstrec,lastrec)
!            read(10,rec=irec) inner_ior4(i1dim:i2dim,j1dim:j2dim)
!            WRITE(66,'(14e15.6)')((inner_ior4(i,j),i=i1,i2),j=j1,j2)
!            !print *,irec,maxval(inner_ior4(i1dim:i2dim,j1dim:j2dim)),i1dim,i2dim,j1dim,j2dim
!            do jx=j1,j2
!            do ix=i1,i2
!               ixx=max(ix-i1,jx-j1)+1
!               mat_iofld(ixx,k)=inner_ior4(ix,jx)
!            end do
!            end do
!         end do
!#if defined(MATLAB)
!         call mxCopyReal8ToPtr(mat_iofld,mxGetPr(pa1),int(imax*kdm,kind=IKIND))
!         iret=matPutVariable(mp, 'vbaclin'//cll, pa1)
!#endif /*MATLAB*/
!
!         ! Get variable
!         do k=1,kdm
!            call get_nest_record('UT ',k,nestingfile,irec,firstrec,lastrec)
!            read(10,rec=irec) inner_ior4(i1dim:i2dim,j1dim:j2dim)
!            tmpfld=inner_ior4
!            call get_nest_record('UB ',0,nestingfile,irec,firstrec,lastrec)
!            read(10,rec=irec) inner_ior4(i1dim:i2dim,j1dim:j2dim)
!            tmpfld=tmpfld+inner_ior4
!            WRITE(66,'(10e15.6)')((tmpfld(i,j),i=i1,i2),j=j1,j2)
!            !print *,irec,maxval(inner_ior4(i1dim:i2dim,j1dim:j2dim)),i1dim,i2dim,j1dim,j2dim
!            do jx=j1,j2
!            do ix=i1,i2
!               ixx=max(ix-i1,jx-j1)+1
!               mat_iofld(ixx,k)=tmpfld(ix,jx)
!            end do
!            end do
!         end do
!#if defined(MATLAB)
!         call mxCopyReal8ToPtr(mat_iofld,mxGetPr(pa1),int(imax*kdm,kind=IKIND))
!         iret=matPutVariable(mp, 'utot'//cll, pa1)
!#endif /*MATLAB*/
!
!         ! Get variable
!         do k=1,kdm
!            call get_nest_record('VT ',k,nestingfile,irec,firstrec,lastrec)
!            read(10,rec=irec) inner_ior4(i1dim:i2dim,j1dim:j2dim)
!            tmpfld=inner_ior4
!            call get_nest_record('VB ',0,nestingfile,irec,firstrec,lastrec)
!            read(10,rec=irec) inner_ior4(i1dim:i2dim,j1dim:j2dim)
!            tmpfld=tmpfld+inner_ior4
!            WRITE(66,'(10e15.6)')((tmpfld(i,j),i=i1,i2),j=j1,j2)
!            !print *,irec,maxval(inner_ior4(i1dim:i2dim,j1dim:j2dim)),i1dim,i2dim,j1dim,j2dim
!            do jx=j1,j2
!            do ix=i1,i2
!               ixx=max(ix-i1,jx-j1)+1
!               mat_iofld(ixx,k)=tmpfld(ix,jx)
!            end do
!            end do
!         end do
!#if defined(MATLAB)
!         call mxCopyReal8ToPtr(mat_iofld,mxGetPr(pa1),int(imax*kdm,kind=IKIND))
!         iret=matPutVariable(mp, 'vtot'//cll, pa1)
!#endif /*MATLAB*/
!
!
!         ! Ice, only matlab for now
!         call get_nest_record('HI ',0,nestingfile,irec,firstrec,lastrec)
!         read(10,rec=irec) inner_ior4(i1dim:i2dim,j1dim:j2dim)
!         tmpfld=inner_ior4
!         do jx=j1,j2
!         do ix=i1,i2
!            ixx=max(ix-i1,jx-j1)+1
!            mat_iofld(ixx,1)=tmpfld(ix,jx)
!         end do
!         end do
!#if defined(MATLAB)
!         pa2=mxCreateDoubleMatrix(int(imax,kind=IKIND),int(1,kind=IKIND),int(0,kind=IKIND))
!         call mxCopyReal8ToPtr(mat_iofld(:,1),mxGetPr(pa2),int(imax,kind=IKIND))
!         iret=matPutVariable(mp, 'hice'//cll, pa2)
!#endif /*MATLAB*/
!
!         ! Ice, only matlab for now
!         call get_nest_record('FI ',0,nestingfile,irec,firstrec,lastrec)
!         read(10,rec=irec) inner_ior4(i1dim:i2dim,j1dim:j2dim)
!         tmpfld=inner_ior4
!         do jx=j1,j2
!         do ix=i1,i2
!            ixx=max(ix-i1,jx-j1)+1
!            mat_iofld(ixx,1)=tmpfld(ix,jx)
!         end do
!         end do
!#if defined(MATLAB)
!         call mxCopyReal8ToPtr(mat_iofld(:,1),mxGetPr(pa2),int(imax,kind=IKIND))
!         iret=matPutVariable(mp, 'fice'//cll, pa2)
!#endif /*MATLAB*/
!
!
!      end do
!
!
!      close(66)
!      close(10)
!
!      print *,' Data dumped to tecplot file '//trim(nestingfile)//'_'//cbnd//'_test.tec'
!
!#if defined(MATLAB)
!      print *,' Data dumped to matlab  file '//trim(nestingfile)//'_'//cbnd//'_test.mat'
!#endif






      end program check_nest
