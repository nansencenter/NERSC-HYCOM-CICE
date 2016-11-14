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
   implicit none
   character(len=8) :: fldextract(1000), char8,dpfield
   integer, parameter :: iversion=4
   integer ivers
   integer :: velrot
   logical :: norot
   logical :: appendfile
   integer, parameter :: maxfld=1000
   character(len=80) :: fnamein, ftype
   type(hycomfile) :: hfile
   character(len=8)  :: ufld(maxfld)
   integer :: nunique, itmp, nfld, i, j, k, kdm
   logical :: option
   real, dimension(:,:), allocatable :: field, dp
   real, dimension(:,:,:), allocatable :: intf

   ! parse arguments can be found at the end of this file
   call parse_arguments()

   ! What file is this? (daily? weekly? restart? pak?)
   ftype=getfiletype(trim(fnamein))

   ! Inits file type
   call initHF(hfile,fnamein,trim(ftype))

   ! determine number of and which fields to be read -- Now uses 
   ! same input files as m2nc
   if(trim(ftype)=="restart") then
      open(21,file='extract.restart',status='old')
   elseif(trim(ftype)=="nersc_daily") then
      open(21,file='extract.daily',status='old')
   elseif(trim(ftype)=="nersc_weekly") then
      open(21,file='extract.weekly',status='old')
   elseif(trim(ftype)=="archv") then
      open(21,file='extract.archv',status='old')
   end if
   read(21,*) ivers  
   if (ivers /= iversion) then
      print *,'wrong version of extract file',ivers,iversion
      call exit(1)
   endif
   read(21,*) ! ignored
   read(21,*) ! ignored
   read(21,*) ! ignored
   read(21,*) ! ignored
   nfld=0
   do i=1,1000
      read(21,*,end=888) char8 ,itmp,itmp,option
      if (option) then
         nfld=nfld+1
         fldextract(nfld)=char8
      end if
   enddo
888 close(21)

   ! TODO - use m_fields_to_plot and vector logic



   ! Initialize IO for .ab files
   CALL XCSPMD  ! -- Requires "regional.grid.b" to be present
   CALL ZAIOST

   ! Get model grid & Depths
   call get_grid()

   ! Read nodes along section -- this assumes "section_intersect" is already run
   call read_section_nodes()

   ! 
   kdm=vDim(hfile)
   allocate(intf(idm,jdm,1:kdm+1))
   allocate(field(idm,jdm))
   allocate(dp(idm,jdm))
   call ncwrite_secdata('lon',plon,0,kdm,var3d=.false.,appendfile=appendfile)
   call ncwrite_secdata('lat',plat,0,kdm,var3d=.false.,appendfile=appendfile)
   call ncwrite_secdata('depth',depths,0,kdm,var3d=.false.,appendfile=appendfile)

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
   ! TODO: Velocity rotation
   ! TODO: metadata
   do j=1,nfld
      if (is3DVar(hfile,fldextract(j),1)) then

         write(6,'(a,a8,x,i1)',advance='no') '3D Var: ',fldextract(j),1

         do k=1,kdm
            write(6,'(a1)',advance='no') '.' ; call flush(6)
            if (k==kdm) write(6,'(i2)',advance='yes') kdm ; call flush(6)

            call HFReadField(hfile,field,idm,jdm,fldextract(j),k,1)

            ! Unmask near-empty layers
            dp=intf(:,:,k+1) - intf(:,:,k)
            where( dp < .1) field=undef

            ! put into netcdf file(s)
            call ncwrite_secdata(fldextract(j),field,k,kdm, &
                 var3d=.true.,vartime=fyear(hfile), appendfile=appendfile, &
                 fillvalue=undef)

         end do
      else



         write(6,'(a)') '2D Var: '//trim(fldextract(j))
         call HFReadField(hfile,field,idm,jdm,fldextract(j),0,1)

         ! put into netcdf file(s) 
         call ncwrite_secdata(fldextract(j),field,0,kdm, &
            vartime=fyear(hfile),appendfile=appendfile, &
            fillvalue=undef)


      end if
   end do





      

   contains


      subroutine parse_arguments()
      implicit none
      integer :: iarg
      character(len=80) :: tmpc
      integer*4, external :: iargc

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
