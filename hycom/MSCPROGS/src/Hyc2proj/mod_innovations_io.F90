module mod_innovations_io

contains
   
   subroutine innovations_write( obslon,obslat, obsdepth,  mod_anomal,mod_mean, obsval,obsvar, obsid, &
   nobs,nrens,filename)
      implicit none

      
      integer,           intent(in) :: nobs,nrens
      character(len=*),  intent(in) :: filename
      real             , intent(in) :: mod_anomal(nobs,nrens),mod_mean(nobs)
      real,              intent(in), dimension(nobs) :: obslon,obslat,obsvar,obsval, obsdepth
      character(len=5),  intent(in), dimension(nobs) :: obsid

      character(len=80) :: frm_string
      character(len=7) :: cens
      integer :: j,iobs,modulus_obs
      logical :: l_analyzed


      modulus_obs = nobs / 8


      ! Length of one record ... 
      ! lon  |  lat  |  nrens  |  id  |  obs  |  obsvar  |  model values
      inquire(iolength=j)                                 &
         obslon(1),obslat(1), obsdepth(1),nrens,obsid(1),obsval(1),obsvar(1), &
         mod_anomal(1,1:nrens)+mod_mean(iobs)

      write (cens,'(i6)') nrens
      frm_string='(3f8.2,i5,a7,2f10.2,'//trim(cens)//'f10.2)'
      !print *,cens,frm_string

      open(10,file=trim(filename),form='unformatted',access='direct',recl=j)
      do iobs=1,nobs
         write(10,rec=iobs)                               &
         obslon(iobs),obslat(iobs),obsdepth(iobs), nrens,obsid(iobs), obsval(iobs),obsvar(iobs), &
         mod_anomal(iobs,1:nrens)+mod_mean(iobs)

         !if (mod(iobs,modulus_obs)==0) then
         !   write(*,trim(frm_string)) &
         !   obslon(iobs),obslat(iobs),obsdepth(iobs),nrens,obsid(iobs), &
         !   obsval(iobs),obsvar(iobs), mod_anomal(iobs,1:nrens)+mod_mean(iobs)
         !end if
      end do
      close(10)
      !print *,cens,frm_string

   end subroutine innovations_write

   





   subroutine innovations_read(lon,lat,depth,modval,obsval,obsvar,obsid,nrobs,nrens,filename)
      implicit none

      real,  dimension(:)  , pointer :: lon, lat,depth
      character(len=5), dimension(:)  , pointer :: obsid
      real, dimension(:)  , pointer :: obsval,obsvar
      real, dimension(:,:), pointer :: modval
      integer,   intent(out) :: nrobs,nrens
      character(len=*), intent(in) :: filename


      character(len=80) :: frm_string
      character(len=7) :: cens
      character(len=5) :: cid
      integer :: j,iobs,ios,itmp, modulus_obs
      logical :: l_analyzed
      real    :: tlon, tlat, tobsval, tobsvar, tdepth
      real, allocatable :: tmpens(:)



      ! First test length of record
      inquire(iolength=j) tlon,tlat,tdepth,nrens
      open(10,file=trim(filename),form='unformatted',access='direct',recl=j,iostat=ios)
      read(10,rec=1) tlon,tlat,tdepth,nrens
      close(10)
      print *,'nrens= ',nrens
      if (ios/=0.or.nrens<1.or.nrens>1e5) then
         print *,'Dubious number of members, or read error in innovations_read'
         print *,'iostat=',ios
         print *,'nrens= ',nrens
         print *,'(If the number of members IS correct, then fix me ....)'
         write(0,'(i2)') 1
         stop '(m_innovations_read.F90)'
      end if
      print *,'filename=',trim(filename)



      ! Second test # records
      allocate(tmpens(nrens))
      inquire(iolength=j) tlon,tlat,tdepth,itmp,cid,tobsval,tobsvar,tmpens
      open(10,file=trim(filename),form='unformatted',access='direct',recl=j)
      ios=0
      iobs=0
      do while (ios==0)
         iobs=iobs+1
         read(10,rec=iobs,iostat=ios)  tlon,tlat,tdepth,itmp,cid,tobsval,tobsvar,tmpens
      end do
      nrobs=iobs-1
      print *,'nrobs= ',nrobs
      if (nrobs<1.or.nrobs>1e7) then
         print *,'Dubious number of obs'
         print *,'nrens= ',nrobs
         print *,'(If the number of obs IS correct, then fix me ....)'
         write (0,'(i2)') 2
         stop '(m_innovations_read.F90)'
      end if

      ! Set this to > nrobs for no output ...
      modulus_obs = nrobs / 8
      

      ! Allocate
      allocate(lon   (nrobs))
      allocate(lat   (nrobs))
      allocate(depth (nrobs))
      allocate(obsval(nrobs))
      allocate(obsvar(nrobs))
      allocate(obsid (nrobs))
      allocate(modval(nrobs,nrens))

      !print *,'cens'
      write (cens,'(i6)') nrens
      !print *,'frm_string'
      frm_string='(3f8.2,i5,a7,2f10.2,'//trim(cens)//'f10.2)'
      print *, 'cens , format_string :',cens,frm_string

      ios=0
      do iobs=1,nrobs
         read(10,rec=iobs,iostat=ios)  &
         lon(iobs),lat(iobs),depth(iobs),itmp,obsid(iobs),&
         obsval(iobs), obsvar(iobs),modval(iobs,:)

         !if (mod(iobs,modulus_obs)==0) then
         !   write(*,trim(frm_string))  &
         !   lon(iobs),lat(iobs),depth(iobs),itmp,obsid(iobs), &
         !   obsval(iobs), obsvar(iobs),modval(iobs,:)
         !end if
      end do
      close(10)
      if (ios/=0) then
         print *,'Read Error ....'
         print *,'ios=',ios
         write(0,'(i2)') 2
         stop '(m_innovations_read.F90)'
      end if


      ! If we came this far.. everything should be ok. Dump 0 to standard error
      write(0,'(i2)') 0




   end subroutine innovations_read

end module mod_innovations_io


