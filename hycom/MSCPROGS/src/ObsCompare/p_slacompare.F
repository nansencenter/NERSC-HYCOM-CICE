      program p_slacompare
      ! Routine to compare modeled sla with observations 
      ! Program takes three arguments  - first is model
      ! field (only DODS link for now), second is SLA
      ! field - also DODS. Third is some form of time 
      ! unit (to be decided)
      !
      ! This program must be linked against nc-dap
      use mod_xc, only: idm,jdm, xcspmd
      use mod_za, only: zaiost
      use mod_sphere_tools
      use mod_grid 
      use m_duacsread
      use mod_confmap
      use mod_year_info
      use mod_hycomfile_io
      use mod_regions
      implicit none

#if defined (IARGC)
      integer*4, external :: iargc
#endif
      !include 'nc-dap.inc'
      include 'netcdf.inc'

      integer, parameter :: nblon=1080 ,     nblat=915 
      real, parameter :: undef2=-1e14 , noasciidata=1e4
      real :: lon(nblon), lat(nblat), sla(nblat, nblon)

      character(len=80) :: hycfile, fbase
      type(hycomfile) :: hfile

      character(len=3) :: creg
      integer, allocatable :: regflag(:,:)
      real, allocatable :: modssh(:,:), modsla(:,:),newsla(:,:)
      logical, allocatable ::regmask(:,:)
      real*8, allocatable :: meanssh(:,:)
      character(len=11) tag7
      logical :: ex,havedata
      integer :: i,j,ipiv,jpiv,jday,ireg
      real :: lat_n,lon_n
      type(year_info) :: rti, rtd
      integer :: stat, numregion
      integer :: ncid, dimid_idm, dimid_jdm, obs_varid, mod_varid, 
     &   mod_msshvarid, mod_sshvarid, maskvarid


      ! Regional stats
      integer, parameter :: maxreg=20, maxcorner=20
      integer :: ios, iosa, iosb, iosc, ierr
      integer :: ncorners
      real, dimension(:), allocatable :: clons, clats
      character(len=80) :: filetype, duacsfile
      integer :: npoints, npoints2
      real :: diffsq, meanobs, meanmod, meanmodonly




c --- Process input args
      if (iargc()/=1 .and. iargc()/=2) then
          print *,'(a)','************ slacmp *************************'
          print *,'(a)','*This routine will compare model SLA with   *'
          print *,'(a)','*SLA from aviso-duacs. If OpenDAP support   *'
          print *,'(a)','*is available, the routine will try to read *'
          print *,'(a)','*SLA data from the opendap address of DUACS *' 
          print *,'(a)','*In this case, the only argument to this    *' 
          print *,'(a)','*routine is a hycom file.                   *' 
          print *,'(a)','*If opendap support is not compiled in, the *' 
          print *,'(a)','*second argument must be a DUACS merged sla *' 
          print *,'(a)','*file retrieved from their web site(merged  *' 
          print *,'(a)','*means a map file, not tracks)              *' 
          print *,'(a)','*                                           *' 
          print *,'(a)','*You will also need a regiondefs.in file    *' 
          print *,'(a)','*for statistics, and a mean ssh file from   *' 
          print *,'(a)','*the model                                  *' 
          print *,'(a)','************ slacmp *************************'
          print *,'Usage:'
          print *,'  slacmp hycom-file [duacs-sla-file]'
          stop
       else
          call getarg(1,hycfile)
          if (iargc()==2) call getarg(2,duacsfile)
       end if

c --- Deduce hycom filetype
      filetype= getfiletype(hycfile)

c --- Initialize the hycomfile module - and retrieve vertical levels
      call initHF(hfile,trim(hycfile),trim(filetype))

c --- get julian day rel 1950-1-1
      jday=datetojulian(hfile%iyear,1,1,1950,1,1) + hfile%iday

c --- Read duacs data from opendap
      if (iargc()==1) then ! this uses opendap data access (flaky)
         print *,'duacsread OpenDAP'
         call duacsread(sla,lon,lat,nblon,nblat,undef2,jday)
      else
         print *,'duacsread NetCDF'
         call duacsread(sla,lon,lat,nblon,nblat,undef2,jday,duacsfile)
      end if
      print *

c --- Get model grid
      call xcspmd()
      call zaiost()
      call get_grid() ;
      call initconfmap(idm,jdm)

c --- Read model ssh field
      allocate(modssh(idm,jdm))
      call HFReadField(hfile,modssh,idm,jdm,'ssh     ',0,1)
      where (depths<.1) modssh=undef2

c --- Read mean model ssh
      allocate(meanssh(idm,jdm))
      if (idm>999 .or. jdm>999) then
         write(tag7,'(i5.5,a,i5.5)')idm,'x',jdm
      else
         write(tag7,'(i3.3,a,i3.3)')idm,'x',jdm
      end if
      inquire(exist=ex, file='meanssh'//trim(tag7)//'.uf')
      if (ex) then
         open(11,file  ='meanssh'//trim(tag7)//'.uf', 
     &           status='old',form='unformatted')
         read (11)meanssh
         close(11)
      else
         print *,'Could not find mean ssh file '
         call exit(1)
      end if
      where (depths<.1) meanssh=undef2

c --- Calculate model sla
      allocate(modsla(idm,jdm))
      modsla=0.
      modsla=modssh-real(meanssh,kind=4)
      where (depths<.1) modsla=undef2

c --- Find model points corr  to sla points.
c --- NB: Only use observation points
      allocate(newsla(idm,jdm))
      newsla=undef2
      do j=1,nblon
      do i=1,nblat
         call oldtonew(lat(i),lon(j),lat_n,lon_n)
         call pivotp(lon_n,lat_n,ipiv,jpiv)
         if ( ipiv>1 .and. ipiv<idm .and. 
     &        jpiv>1 .and. jpiv<jdm ) then
            if (depths(ipiv,jpiv)>1. ) then
               newsla(ipiv,jpiv)=sla(i,j)   
            end if
         end if
      end do
      end do


c --- Quick dump of data points to netcdf file
      stat = NF_CREATE('slacmp.nc', NF_CLOBBER , ncid)
      if (stat /= nf_noerr) then
         print *,nf_strerror(stat)
         call exit(3)
      end if
      stat = NF_DEF_DIM (NCID, 'idm', idm, dimid_idm)
      stat = NF_DEF_DIM (NCID, 'jdm', jdm, dimid_jdm)

      stat = NF_DEF_VAR(NCID, 'slaobs', NF_REAL, 2,
     &   (/dimid_idm,dimid_jdm/), obs_varid)

      stat = NF_DEF_VAR(NCID, 'slamod', NF_REAL, 2,
     &   (/dimid_idm,dimid_jdm/), mod_varid)

      stat = NF_DEF_VAR(NCID, 'modmeanssh', NF_REAL, 2,
     &   (/dimid_idm,dimid_jdm/), mod_msshvarid)

      stat = NF_DEF_VAR(NCID, 'modssh', NF_REAL, 2,
     &   (/dimid_idm,dimid_jdm/), mod_sshvarid)

      stat = NF_PUT_ATT_REAL  (NCID, obs_VARID, '_FillValue',
     &   NF_REAL,1, real(undef2,kind=4) )
      stat = NF_PUT_ATT_REAL  (NCID, mod_VARID, '_FillValue',
     &   NF_REAL,1, real(undef2,kind=4) )
      stat = NF_PUT_ATT_REAL  (NCID, mod_msshVARID, '_FillValue',
     &   NF_REAL,1, real(undef2,kind=4) )
      stat = NF_PUT_ATT_REAL  (NCID, mod_sshVARID, '_FillValue',
     &   NF_REAL,1, real(undef2,kind=4) )
      stat = NF_PUT_ATT_REAL  (NCID, obs_VARID, 'missing_value',
     &   NF_REAL,1, real(undef2,kind=4) )
      stat = NF_PUT_ATT_REAL  (NCID, mod_VARID, 'missing_value',
     &   NF_REAL,1, real(undef2,kind=4) )
      stat = NF_PUT_ATT_REAL  (NCID, mod_msshVARID, 'missing_value',
     &   NF_REAL,1, real(undef2,kind=4) )
      stat = NF_PUT_ATT_REAL  (NCID, mod_sshVARID, 'missing_value',
     &   NF_REAL,1, real(undef2,kind=4) )

      STAT = NF_ENDDEF(NCID)


      stat = NF_PUT_VARA_REAL(NCID,mod_VARID,(/1,1/),(/idm,jdm/),
     &   real(modsla,kind=4))
      if (stat /= nf_noerr) then
         print *,nf_strerror(stat)
         call exit(3)
      end if

      stat = NF_PUT_VARA_REAL(NCID,mod_msshVARID,(/1,1/),(/idm,jdm/),
     &   real(meanssh,kind=4))
      if (stat /= nf_noerr) then
         print *,nf_strerror(stat)
         call exit(3)
      end if

      stat = NF_PUT_VARA_REAL(NCID,obs_VARID,(/1,1/),(/idm,jdm/),
     &   real(newsla,kind=4))
      if (stat /= nf_noerr) then
         print *,nf_strerror(stat)
         call exit(3)
      end if


      stat = NF_PUT_VARA_REAL(NCID,mod_sshVARID,(/1,1/),(/idm,jdm/),
     &   real(modssh,kind=4))
      if (stat /= nf_noerr) then
         print *,nf_strerror(stat)
         call exit(3)
      end if






!!!!!!! Format for regional stats file (regiondefs.in):
!ncorners regname
!clons
!clats
!!!!!!!
! Example: 
!4 NORTHATLANTIC
!
!-80 20 20 -80
!40  40 60  60
!!!!!!!
! You can have several region defs in one file, if they follow eachother



c --- Initialize regions for RMS calc
      havedata=.not.all(modssh==undef2)
      call read_regions()
      numregion=getnregions()
      allocate(regmask(idm,jdm))
      allocate(regflag(idm,jdm))

c --- Loop over regions and calculate statistics
      do ireg=1,numregion
         print *
         print '(a,i3,a)','Region ',ireg,' '//trim(getname(ireg))
         npoints=0
         npoints2=0
         diffsq=0.
         meanobs=0.
         meanmod=0.
         meanmodonly=0.

c ---    Get mask for region
         call getmask(ireg,regmask,plon,plat,idm,jdm)


         do j=1,jdm
         do i=1,idm
            if (regmask(i,j) .and. depths(i,j)>1. .and. 
     &         newsla(i,j)/=undef2) then
               npoints=npoints+1
               meanobs=meanobs+newsla(i,j)
               meanmod=meanmod+modsla(i,j)*100.
               diffsq=diffsq+(newsla(i,j)-modsla(i,j)*100.)**2
            end if

            if (regmask(i,j) .and. depths(i,j)>1 ) then
               npoints2=npoints2+1
               meanmodonly=meanmodonly+modsla(i,j)*100.
            end if

c ---       Binary mask (ex. 8 regions = binary 11111111)
c ---       Used for fiagnostics
            if (regmask(i,j)) then
               regflag(i,j)=regflag(i,j)+2**(ireg-1)
            end if
         end do
         end do

         if (npoints2>0) then
            meanmodonly=meanmodonly/npoints2
         end if

c ---    Put RMS/bias data into ascii files
         if (.not.havedata) then
            open(11,file='SLA.RMS.'//trim(getname(ireg))//'.asc',
     &              position='append',status='unknown')
            write(11,'(5f14.5,i8)')  
     &          hfile%fyear,
     &         noasciidata,
     &         noasciidata,
     &         noasciidata,
     &                  0
            close(11)
         elseif (npoints>0) then
            print *,'# points ',npoints
            print *,'mean obs ',meanobs/npoints
            print *,'mean mod ',meanmod/npoints
            print *,'bias     ',(meanmod-meanobs)/npoints
            print *,'RMS      ',sqrt(diffsq/npoints)
            open(11,file='SLA.RMS.'//trim(getname(ireg))//'.asc',
     &              position='append', status='unknown')
            write(11,'(5f14.5,i8)')
     &         hfile%fyear,
     &         hfile%fyear,
     &         meanmodonly,              
     &         (meanmod-meanobs)/npoints,
     &         sqrt(diffsq/npoints), 
     &         npoints
            close(11)
         end if
      end do

c --- put mask in netcdf file for checks
      if (ireg>0) then
         STAT = NF_REDEF(NCID)
         stat = NF_DEF_VAR(NCID, 'regionmask',NF_INT,
     &      2,(/dimid_idm,dimid_jdm/), maskvarid)
         if (stat /= nf_noerr) then
            print *,nf_strerror(stat)
            call exit(3)
         end if
         STAT = NF_ENDDEF(NCID)

         where(depths<.1) regflag=-2
         stat = NF_PUT_VARA_INT(NCID,maskvarid,
     &             (/1,1/),(/idm,jdm/),regflag)
         if (stat /= nf_noerr) then
            print *,nf_strerror(stat)
            call exit(3)
         end if
      end if
      stat = NF_CLOSE(ncid)

c --- Final info
      print *
      print '(a)','Output in '
      do ireg=1,numregion
         print *,'SLA.RMS.'//trim(getname(ireg))//'.asc'
      end do
      print *,'Fields in slacmp.nc'



      end program

