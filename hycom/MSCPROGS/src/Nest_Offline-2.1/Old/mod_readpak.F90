module mod_readpak

contains

subroutine get_pakfield(filename,fldname,klevel,twod,idm,jdm,undef)
   implicit none
   character(len=*), intent(in) :: filename, fldname
   integer, intent(in) :: klevel,idm,jdm
   real,    intent(out):: twod(idm,jdm)
   real, intent(in) :: undef


   real :: tmpwod(idm,jdm)
   character(len=2) :: util(idm*jdm+14)
   character(len=80)  fieldhead
   character(len=5) :: pakfldname,recfldname
   logical :: found
   integer ::n,ios,j,dataunit,reclayer

   dataunit=33


   !if (klevel==1) then
   !   write(*,'(a,a5)',advance='no') 'looking for variable ',fldname
   !end if

   ! open
   inquire(iolength=j)fieldhead,util
   open(dataunit,file=trim(filename),form='unformatted',access='direct', &
        action='read',recl=j)

   ! Convert fieldname to pak fieldname
   call getpakname(fldname,pakfldname)

   found=.false.
   ios=0


   ! Find field header corresponding to fld
   print *,fldname,pakfldname,klevel
   tmpwod=0.
   n=1
   do while(ios==0.and..not.found)
      read(dataunit,rec=n,iostat=ios)fieldhead,util
      read(fieldhead,'(t4,a5,i2)') recfldname,reclayer
      !print *,recfldname,reclayer

      if (trim(adjustl(recfldname))==trim(pakfldname) .and. reclayer==klevel) then
         print *,'|---->found ok '
         call unpakk(twod,idm,idm,jdm,util,idm*jdm+14)
         found=.true.
      endif
      n=n+1
   enddo

   if (.not.found) then
   !   write(*,*)
   !   write(*,'(a,a5,a,i2,a)')'DID NOT FIND',&
   !         hdr(n)%fieldname,'(',hdr(n)%layer,')'
   !   write(*,'(a,a5,a,i2,a)')'DID NOT FIND',&
   !         fld%fextract,'(',k,')'
   !   stop
   !else if (k==1) then
   !   write(*,'(a)') '....... ok'
      twod=undef
   endif

   ! Interpolate velocity fields to p-grid
   !if (hdr(n)%fieldname(1:1)=='U') then
   !   tmpwod=0.
   !   tmpwod(1:idm-1,:) = .5*(twod(1:idm-1,:)+twod(2:idm,:))
   !   twod = tmpwod
   !else if (hdr(n)%fieldname(1:1)=='V') then
   !   tmpwod=0.
   !   tmpwod(:,1:jdm-1) = .5*(twod(:,1:jdm-1)+twod(:,2:jdm))
   !   twod = tmpwod
   !end if

   close(dataunit)

   end subroutine get_pakfield




subroutine read_pakheader(fname,nx,ny,nz,runcode,plot_time,c2d,c3d,n2d,n3d,maxc)
   implicit none
   integer, intent(out) :: nx
   integer, intent(out) :: ny
   integer, intent(out) :: nz
   integer, intent(out) :: n2d
   integer, intent(out) :: n3d
   integer, intent(in) :: maxc
   character(len=5), intent(out), dimension(maxc)  ::  c2d,c3d
   character(len=3), intent(out) :: runcode
   !character(len=24), intent(out) :: plot_time
   character(len=*), intent(out) :: plot_time
   character(len=*), intent(in) :: fname

   character(len=50) line
   character(len=11) char11
   character(len=7) char7
   character(len=5) char5
   integer nstep,i


! Read header from headerfile
   print '(a)','#############################################'
   print *,'Info from plot header file :'
   open(10,file=trim(fname),action='read')
      read(10,'(a50)')line 
      if (line(1:3) == '1.1') then
         print *,'header version 1.1'
         read(10,'(t6,a3)')runcode
         !read(10,'(a24)')plot_time
         read(10,'(a)')plot_time
      else
         plot_time(:)=' '
         plot_time(1:12)='time unknown'
         runcode(1:3)='RUN'
         print *,'old header file'
         read(10,'(A7,I8)') char7,nstep
         read(10,*)
      endif
      read(10,'(A11,I2)')char11,n2d
      print *,'2D fields in dump file = ',n2d
      do i=1,n2d
         read(10,'(a5)')char5
         c2d(i)=char5
         write(*,'(a,a)',advance='no')char5,','
      enddo
      write(*,*)
      read(10,'(A11,I2)')char11,n3d
      print *,'3D fields in dump file = ',n3d
      do i=1,n3d
         read(10,'(a5)')char5
         c3d(i)=char5
         write(*,'(a,a)',advance='no')char5,','
      enddo
      write(*,*)
      read(10,'(t11,I4,t18,i4,t25,i4)')nx,ny,nz
   close(10)
   print '(a,3I4)','Grid dimensions are : ',nx,ny,nz
   print '(a,a3)', 'Run code is         : ',runcode
   print '(a,a24)','Plot time is        : ',plot_time
   print '(a)','#############################################'
end subroutine read_pakheader


   ! Conversion table for pak field names
   subroutine getpakname(fldname,pakfldname)
   implicit none
   character(len=*), intent(in)  :: fldname
   character(len=*), intent(out) :: pakfldname

   select case (trim(fldname)) 
      case('temp')
         pakfldname='TEM'
      case('saln')
         pakfldname='SAL'
      case('utot')
         pakfldname='UT'
      case('vtot')
         pakfldname='VT'
      case('pres')
         pakfldname='DP'
      case('ubavg')
         pakfldname='UBAVG'
      case('vbavg')
         pakfldname='VBAVG'
      case default
         pakfldname=''
   end select
   end subroutine


end module mod_readpak
