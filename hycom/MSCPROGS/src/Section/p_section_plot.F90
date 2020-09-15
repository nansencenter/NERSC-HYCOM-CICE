! Program to convert from various hycom - type files
! to netcdf files with data interploated along a section
! KAL: Added options to read other than pak files
!

program section_plot
   use mod_xc
   use mod_za
   use mod_sections
   use mod_grid
   use mod_hycomfile_io
   use m_fields_to_plot
   implicit none
   character(len=8) :: char8
   integer ivers
   integer :: velrot
   logical :: norot
   logical :: appendfile
   integer, parameter :: maxfld=1000
   character(len=80) :: fnamein, ftype
   type(hycomfile) :: hfile
   character(len=20)  :: uname,vname
   integer :: nunique, itmp,  i, j, k, kdm
   logical :: option
   real, dimension(:,:), allocatable :: field, dp, field2
   real, dimension(:,:,:), allocatable :: intf
   type(fields) :: fld(1000)
   integer :: nfld

   ! parse arguments can be found at the end of this file
   call parse_arguments()

   write(6,'(a)')
   ! What file is this? (daily? weekly? restart? pak?)
   ftype=getfiletype(trim(fnamein))
   print *, ftype
   
   write(6,'(a)')
   ! Inits file type
   call initHF(hfile,fnamein,trim(ftype))

   ! Gets fields to plot, assuming ftype is same for all files (!)
   call fields_to_plot(fld,nfld,hfile)

   ! Initialize IO for .ab files
   CALL XCSPMD  ! -- Requires "regional.grid.b" to be present
   CALL ZAIOST

   ! Get model grid & Depths
   call get_grid()
   
   write(6,'(a)')
   
   ! Read nodes along section -- this assumes "section_intersect" is already run
   call read_section_nodes()

   ! 
   kdm=vDim(hfile)
   allocate(intf(idm,jdm,1:kdm+1))
   allocate(field (idm,jdm))
   allocate(field2(idm,jdm)) ! For vectors
   allocate(dp(idm,jdm))
   call ncwrite_secdata('lon',plon,0,kdm,var3d=.false.,appendfile=appendfile)
   call ncwrite_secdata('lat',plat,0,kdm,var3d=.false.,appendfile=appendfile)
   call ncwrite_secdata('depth',depths,0,kdm,var3d=.false.,appendfile=appendfile)
   !call ncwrite_secdata('dist',ndedist,0,kdm,var3d=.false.,appendfile=appendfile)

   ! We always process the layer interfaces
   dp=0.0
   intf =0.0
   do k=1,kdm
      call HFReadDPField_m(hfile,dp,idm,jdm,k,1) ! Returns thickness in meters
      intf(:,:,k+1)=intf(:,:,k)+dp;
      !write(6,'(a)') '3D Var: intf'
      call ncwrite_secdata('intf_lower',intf(:,:,k+1),k,kdm, &
         vartime=fyear(hfile),var3d=.true.,appendfile=appendfile, &
         comment='Lower interface of layer')
   end do


   !Go through fields -- find  -- interpolate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! Scalar case
   do j=1,nfld
      if (is3DVar(hfile,fld(j)%fextract,1) .and. fld(j)%option) then

         if  (.not. fld(j)%vecflag) then
            write(6,'(a,2a8,x,i1)',advance='no') '3D Scalar Var: ',fld(j)%fextract,''
         else 
            write(6,'(a,2a8,x,i1)',advance='no') '3D Vector Var: ',fld(j)%fextract,fld(j+1)%fextract
         end if
         do k=1,kdm
            write(6,'(a1)',advance='no') '.' ; call flush(6)
            if (k==kdm) write(6,'(i2)',advance='yes') kdm ; call flush(6)

            dp=intf(:,:,k+1) - intf(:,:,k)

            ! Unmask near-empty layers
            call HFReadField(hfile,field,idm,jdm,fld(j)%fextract,k,1)

            ! Read 2nd component for vectors
            if (fld(j)%vecflag) call HFReadField(hfile,field2,idm,jdm,fld(j+1)%fextract,k,1)

            ! Rotate to east/north components for vectors
            if (fld(j)%vecflag) then
               if (velrot==1) then
                  uname='east_'//trim(fld(j)%fextract)
                  vname='north_'//trim(fld(j+1)%fextract)
                  call rotate(field,field2,plat,plon,idm,jdm,'m2l')
               else if (velrot==2) then
                  print *,'Normal rotation not implemented yet'
                  stop
               else if (velrot==0) then
                  uname=trim(fld(j)%fextract)
                  vname=trim(fld(j+1)%fextract)
               else
                  print *,'Unknown rotation flag',velrot
                  stop
               end if
            end if

            ! mask
            where( dp < .1) field=undef
            where( dp < .1) field2=undef

            ! put into netcdf file(s)
            if (.not.fld(j)%vecflag) then ! scalar case
               call ncwrite_secdata(fld(j)%fextract,field,k,kdm, &
                    var3d=.true.,vartime=fyear(hfile), appendfile=appendfile, &
                    fillvalue=undef)

            else ! vector case
               call ncwrite_secdata(trim(uname),field ,k,kdm, &
                    var3d=.true.,vartime=fyear(hfile), appendfile=appendfile, &
                    fillvalue=undef)
               call ncwrite_secdata(trim(vname),field2,k,kdm, &
                    var3d=.true.,vartime=fyear(hfile), appendfile=appendfile, &
                    fillvalue=undef)
            end if
         end do

      elseif (.not.is3DVar(hfile,fld(j)%fextract,1) .and. fld(j)%option) then



         if  (.not. fld(j)%vecflag) then
            write(6,'(a,2a8,x,i1)',advance='no') '2D Scalar Var: ',fld(j)%fextract,''
         else 
            write(6,'(a,2a8,x,i1)',advance='no') '2D Vector Var: ',fld(j)%fextract,fld(j+1)%fextract
         end if

         call HFReadField(hfile,field,idm,jdm,fld(j)%fextract,0,1)
         ! Read 2nd component for vectors
         if (fld(j)%vecflag) call HFReadField(hfile,field2,idm,jdm,fld(j+1)%fextract,0,1)

         ! Rotate to east/north components for vectors
         if (fld(j)%vecflag) then
            if (velrot==1) then
               uname='east_'//trim(fld(j)%fextract)
               vname='north_'//trim(fld(j+1)%fextract)
               call rotate(field,field2,plat,plon,idm,jdm,'m2l')
            else if (velrot==2) then
               print *,'Normal rotation not implemented yet'
               stop
            else if (velrot==0) then
               uname=trim(fld(j)%fextract)
               vname=trim(fld(j+1)%fextract)
            else
               print *,'Unknown rotation flag',velrot
               stop
            end if
         end if

         ! Mask
         where( depths < .1) field=undef
         where( depths < .1) field2=undef

         ! put into netcdf file(s) 
         if (.not.fld(j)%vecflag) then ! scalar case
            call ncwrite_secdata(fld(j)%fextract,field,0,kdm, &
               vartime=fyear(hfile),appendfile=appendfile, &
               fillvalue=undef)

            else ! vector case
            call ncwrite_secdata(uname,field,0,kdm, &
               vartime=fyear(hfile),appendfile=appendfile, &
               fillvalue=undef)
            call ncwrite_secdata(vname,field2,0,kdm, &
                 vartime=fyear(hfile), appendfile=appendfile, &
                 fillvalue=undef)
         end if

      end if
   end do





      

   contains


      subroutine parse_arguments()
      implicit none
      integer :: iarg
      character(len=80) :: tmpc
#if defined(IARGC)
      integer*4, external :: iargc
#endif

      appendfile         =.false.
      velrot=0
      if (iargc()==1) then
         call getarg(1,fnamein)
      elseif (iargc()>1 .and. iargc()<=3) then
         call getarg(1,fnamein)
         do iarg=2,iargc()
            call getarg(iarg,tmpc)
            if (trim(tmpc)=='-append') then
               appendfile=.true.
            elseif (tmpc(1:4)=='-rot') then
               if     (trim(tmpc)=='-rotll') then
                  velrot=1
               elseif (trim(tmpc)=='-rotnormal') then
                  velrot=2
                  print *,'Normal rotation not implemented yet'
                  stop
               else
                  print *,'Allowed rotate arguments are -rotnormal or -rotll'
                  print *, '(section_plot)'
                  call exit(1)
               end if
            else
               print *,'Allowed arguments are -append and -rotnormal or -rotlonlat'
               print *, '(section_plot)'
               call exit(1)
            end if
            end do
      else
         print *,'One or two input arguments only'
         print *, '(section_transport)'
         call exit(1)
      end if
      print *,'Append to file is    ',appendfile
      print *,'Rotate flag is       ',velrot
   end subroutine

   end program section_plot
