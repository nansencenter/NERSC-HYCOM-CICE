module mod_transport

integer, parameter :: maxntrans=300
!We set the max number of criteria to 6 
integer, parameter :: maxncrit=6
integer, save      :: ntrans
character(len=*), parameter :: transfile='transport.in'


#if defined (SCALAR_TRANS) 
character(len=*), parameter :: scalar_transfile='scalartransport.in'
integer, save :: num_scalar_trans               ! # Fields to calc scalar transport for
real, allocatable, dimension(:,:,:) :: scalar   ! Keeps scalar fields
real, allocatable, dimension(:,:)   :: &
   scalar_trans                                 ! Scalar transport values
real, allocatable, dimension(:)     :: &
   scalar_trans_offset, &                       ! Offset, or "zero" value for transport
   scalar_trans_factor                           ! Scalar factor, usually 1
character(len=8), allocatable, dimension(:) :: &
   scalar_trans_name                            ! Name of variable
#endif

character(len=20),save :: transport_name(maxntrans)
integer          ,save :: transport_secnum(maxntrans)
character(len=12),save :: transport_rangeid(maxntrans,maxncrit)
character(len= 1),save :: transport_directionid(maxntrans)
real             ,save :: transport_lowrange(maxntrans,maxncrit)
real             ,save :: transport_uprange (maxntrans,maxncrit)
!Fanf: add an option to use several criteria, nb_crit is the number of criteria
integer          ,save :: transport_nbcrit(maxntrans)

contains


subroutine transport_specification()
use mod_sections
implicit none

   logical :: ex
   integer :: ios,isec,icrit
   character(len=20) :: sectionid
   character(len=12) :: rangeid
   character(len= 1) :: directionid
   real :: lowlim,uplim
   character(len=20) :: transportid
   character(len=80) :: a80
   logical :: eof


   ! Read an infile, specifying density/depth/salinity/temperature intervals
   ! To calculate transport for
   ! Format: '(a20,1x,a10,1x,a1,1x,2f10.2)'
   ! where the fields are: 
   ! 1) Section ID (from sections.in) -- a20
   ! 2) Range ID (a10)  with possibilities
   !    'DEPTH'       -- Specify a depth range
   !    'SALINITY'    -- Specify a salinity range
   !    'TEMPERATURE' -- Specify a temperature range
   !    'DENSITY'     -- Specify a density range
   !    'NONE'        -- No ranges -- use full water column
   ! 3) DIRECTION (a1)
   !    '+' Positive transport across section (in range)
   !    '-' Negative transport across section (in range)
   !    ' ' Net positive transport through section
   ! 4) Two (2f10.2) numbers specifying transport range (see 2)

   inquire(exist=ex,file=transfile)
   if (.not.ex) then
      print *,'Transport file "'//transfile//'" does not exist'
      print *,'Sample input below (format does matter !)'
      write(*,100) 'TRANSPORTID                    ', &
          'SECTIONNAME           ','RANGESPEC      ','+',0.0,10.0
      call exit(1)
   end if

   open(10,file=transfile)
   ios=0
   ntrans=0
   do while (ios==0) 

      eof=.true.
      read(10,100,iostat=ios,end=677) transportid,sectionid,rangeid,directionid,lowlim,uplim
      eof=.false.
677   continue
!Fanf All this part is made for Multiple criteria
!If so then we start a line in transport.in with MULTIPLE
!A new line for each criteria and end up the list with END
      if (trim(rangeID)=='MULTIPLE') then
         icrit=1
         ntrans=ntrans+1
         eof=.true.
         read(10,100,iostat=ios,end=677) transportid,sectionid,rangeid,directionid,lowlim,uplim
         eof=.false.
         do while (trim(rangeID).ne.'END' .and. icrit.le.maxncrit) 
            if (ios==0) then
               write(*,100) transportid,sectionid,rangeid,directionid,lowlim,uplim
               ! Check section ID against available sections
               if (icrit==1) then
                  transport_secnum(ntrans)=-1
                  do isec=1,nsec
                     if (trim(sectionid)==trim(secname(isec))) then
                        transport_secnum(ntrans)=isec
                     end if
                  end do
                  if (transport_secnum(ntrans)==-1) then
                     print *,'No Section corresponding to ID '//trim(sectionid)
                     call exit(1)
                  end if
                  ! Check Direction ID
                  if (trim(directionid)=='+' .or. trim(directionid)=='-' .or.  &
                      trim(directionid)=='n') then 
                     transport_directionid(ntrans)=trim(directionid)
                  else
                     print *,'No direction specification corresponding to '//trim(directionid)
                     call exit(1)
                  end if
               elseif ((trim(sectionid).ne.trim(secname(transport_secnum(ntrans))).or.&
                        trim(transport_directionid(ntrans)).ne.trim(directionid))) then
                     print *,'When MULTIPLE criteria initiated, it should be a  &
                     END before moving to a different SECTIONID or change &
                     transport direction'
                     print *,'icrit is ',icrit
                     print *,' sec name is ',trim(secname(transport_secnum(ntrans))),trim(sectionid)
                     print *,'Trans dir',trim(transport_directionid(ntrans)),trim(directionid)
                     call exit(1)
               end if
               ! Check range ID 
               if (trim(rangeID)=='DEPTH'      .or.trim(rangeID)=='SALINITY'.or. &
                   trim(rangeID)=='TEMPERATURE'.or.trim(rangeID)=='DENSITY' .or. &
                   trim(rangeID)=='NONE' ) then
                  transport_rangeid(ntrans,icrit)=trim(rangeID)
               else
                  print *,'No range specification corresponding to |'//trim(rangeid)//'|'
                  call exit(1)
               end if
               ! Get ranges
               transport_name(ntrans)=transportid
               transport_lowrange(ntrans,icrit)=lowlim
               transport_uprange (ntrans,icrit)=uplim
               transport_nbcrit(ntrans)=icrit
            elseif (ios/=0 .and. .not.eof) then
               print *,'non-eof error: iostat is ',ios
               backspace(10)
               read(10,'(a80)') a80
               print *,'Error in section specification on line ',ntrans
               print *,'Format is in 1st line(| denotes first letter), input line in 2nd line.'
               print *,'Also make sure number positions (decimals) match!'
               write(*,100) '|                   ', &
                            '|                   ', &
                            '|           ',         &
                            '|',                    &
                            0.0,0.0
               write(*,'(a)') a80
               print *, '(transport_specification)'
               call exit(1)
            end if
            eof=.true.
            read(10,100,iostat=ios,end=677) transportid,sectionid,rangeid,directionid,lowlim,uplim
            eof=.false.
            icrit=icrit+1
         end do
      else !SINGLE CRITERIA
         if (ios==0) then
            write(*,100) transportid,sectionid,rangeid,directionid,lowlim,uplim

            ntrans=ntrans+1

            ! Check section ID against available sections
            transport_secnum(ntrans)=-1
            do isec=1,nsec
               if (trim(sectionid)==trim(secname(isec))) then
                  transport_secnum(ntrans)=isec
               end if
            end do
            if (transport_secnum(ntrans)==-1) then
               print *,'No Section corresponding to ID '//trim(sectionid)
               call exit(1)
            end if

            ! Check range ID 
            if (trim(rangeID)=='DEPTH'      .or.trim(rangeID)=='SALINITY'.or. &
                trim(rangeID)=='TEMPERATURE'.or.trim(rangeID)=='DENSITY' .or. &
                trim(rangeID)=='NONE' ) then

               transport_rangeid(ntrans,1)=trim(rangeID)
            else
               print *,'No range specification corresponding to |'//trim(rangeid)//'|'
               call exit(1)
            end if

            ! Check Direction ID
            if (trim(directionid)=='+' .or. trim(directionid)=='-' .or.  &
                trim(directionid)=='n') then 
               transport_directionid(ntrans)=trim(directionid)
            else
               print *,'No direction specification corresponding to '//trim(directionid)
               call exit(1)
            end if

            ! Get ranges
            transport_name(ntrans)=transportid
            transport_lowrange(ntrans,1)=lowlim
            transport_uprange (ntrans,1)=uplim
            transport_nbcrit(ntrans)=1
         elseif (ios/=0 .and. .not.eof) then
            print *,'non-eof error: iostat is ',ios
            backspace(10)
            read(10,'(a80)') a80
            print *,'Error in section specification on line ',ntrans
            print *,'Format is in 1st line(| denotes first letter), input line in 2nd line.'
            print *,'Also make sure number positions (decimals) match!'
            write(*,100) '|                   ', &
                         '|                   ', &
                         '|           ',         &
                         '|',                    &
                         0.0,0.0
            write(*,'(a)') a80
            print *, '(transport_specification)'
            call exit(1)
         end if
      end if  !MULTIPLE/SINGLE CRITERION
   end do
   close (10)


100 format(a20,1x,a20,1x,a12,1x,a1,1x,2f10.2)

end subroutine transport_specification





#if defined (SCALAR_TRANS) 
subroutine scalar_transport_specification(idm,jdm)
   implicit none
   integer, intent(in) :: idm,jdm
   logical :: ex
   character(len=8) :: stname
   real :: stoffset, stfac
   integer :: ios,ist

   ios=0

   inquire(file=scalar_transfile,exist=ex)
   num_scalar_trans=0
   if (ex) then
      open(10,file=scalar_transfile,form='formatted')
      do while (ios==0) 
         read(10,*,iostat=ios) stname, stoffset, stfac
         if (ios==0) num_scalar_trans=num_scalar_trans + 1
      end do
      print *,'Num scalar transports : ',num_scalar_trans
      allocate(scalar(idm,jdm,num_scalar_trans))
      allocate(scalar_trans(ntrans,num_scalar_trans))
      allocate(scalar_trans_name(num_scalar_trans))
      allocate(scalar_trans_offset(num_scalar_trans))
      allocate(scalar_trans_factor(num_scalar_trans))

      scalar_trans=0.

      if (num_scalar_trans>0) then
         rewind(10)
         do ist=1,num_scalar_trans
            read(10,*,iostat=ios) scalar_trans_name(ist), &
               scalar_trans_offset(ist), scalar_trans_factor(ist)
            print '(a,i3,a,a8,a,e14.2,a,e14.2)', ' Scalar transport ',ist,' variable:', &
               scalar_trans_name(ist), ' offset: ', scalar_trans_offset(ist),  &
               ' scale factor:', scalar_trans_factor(ist) 
            if (ios/=0)  then
               print *,'Error readig scalar transport '
               stop '(mod_transport:scalar_transport_specification)'
            end if
         end do
      end if
      close(10)
   end if
   !num_scalar_trans=2
   !allocate(scalar(idm,jdm,num_scalar_trans))
   !allocate(scalar_trans(ntrans,num_scalar_trans))
   !allocate(scalar_trans_name(num_scalar_trans))
   !allocate(scalar_trans_offset(num_scalar_trans))
   !allocate(scalar_trans_sfac(num_scalar_trans))
   !scalar_trans_name(1)='temp'
   !scalar_trans_name(2)='saln'
   !scalar_trans_offset=0.
   !scalar_trans_sfac=0.
end subroutine
#endif











   ! Simplest transport calculation - actually belongs more to the
   ! sections module
   subroutine transport(hfile,appendfile)
      use mod_grid
      use mod_sections
      use mod_hycomfile_io
      implicit none
      type(hycomfile), intent(in) :: hfile
      logical, intent(in), optional :: appendfile
      real   , parameter :: epsiloon=1e-7
      integer :: n,i,j,k,isec, ipiv, jpiv, flagu, flagv, ipnt, ib
      real, dimension(idm,jdm) :: fice, hice,uice,vice,volu,volv, &
         ubavg,vbavg
      real, dimension(max_sec) :: voltr, voltrp, voltrm,  &
         ivoltr, ivoltrm, ivoltrp
      real facu, facv, rtime, rtime2
      character(len=5) :: char5
      character(len=3) :: css
      character(len=80) :: fname
      character(len=40) :: timeunits
      logical :: ex,vapp

      ! set to true to use old text files
      logical, parameter :: textfileflag=.false. 

      rtime=fyear(hfile)

      vapp=.false.
      if (present(appendfile)) vapp=appendfile

      ! Get barotropic velocities
      call HFReaduvbaro(hfile,ubavg,vbavg,idm,jdm,1)

      ! Get ice velocities
         call HFReadField(hfile,uice,idm,jdm,'uice    ',0,1)
         call HFReadField(hfile,vice,idm,jdm,'vice    ',0,1)
      if (trim(hfile%ftype)=='nersc_daily') then
         call HFReadField(hfile,hice,idm,jdm,'hice    ',0,1)
         call HFReadField(hfile,fice,idm,jdm,'fice    ',0,1)
      else if (trim(hfile%ftype)=='nersc_weekly') then
         call HFReadField(hfile,hice,idm,jdm,'hicem   ',0,1)
         call HFReadField(hfile,fice,idm,jdm,'ficem   ',0,1)
      endif
      do j=2,jdm
      do i=1,idm
         if (periodic) then
            ib = mod(idm+i-2,idm)+1
         else
            ib = max(i-1,1)
         end if
         volu(i,j)=0.5*(hice(i,j)*fice(i,j) + hice(i-1,j)*fice(i-1,j))
         volv(i,j)=0.5*(hice(i,j)*fice(i,j) + hice(i,j-1)*fice(i,j-1))
      end do
      end do

      if (textfileflag) then 
         open(345,file='section_transport.filelist',status='replace')
      end if

      voltr =0.
      voltrp=0.
      voltrm=0.
      do isec=1,nsec

         do ipnt=1,sdm(isec)
            ipiv=ndeipiv(ipnt,isec)
            jpiv=ndejpiv(ipnt,isec)
            flagu=ndeflagu(ipnt,isec)
            flagv=ndeflagv(ipnt,isec)
!            print *,isec, ipnt,ipiv,jpiv,flagu,flagv

            ! Total transport
            voltr(isec) = voltr(isec)                                     +  &
                 depthu(ipiv,jpiv)*flagu*ubavg(ipiv,jpiv)*scuy(ipiv,jpiv) +  &
                 depthv(ipiv,jpiv)*flagv*vbavg(ipiv,jpiv)*scvx(ipiv,jpiv) 
            ivoltr(isec) = ivoltr(isec)                                  +  &
                 volu(ipiv,jpiv)*flagu*uice(ipiv,jpiv)*scuy(ipiv,jpiv) +  &
                 volv(ipiv,jpiv)*flagv*vice(ipiv,jpiv)*scvx(ipiv,jpiv) 

            ! Positive transport - contributes if velocity X direction flag > 0
            facu=step(ubavg(ipiv,jpiv)*flagu,0.0)
            facv=step(vbavg(ipiv,jpiv)*flagv,0.0)

            voltrp(isec) = voltrp(isec)                                          +  &
                 depthu(ipiv,jpiv)*flagu*facu*ubavg(ipiv,jpiv)*scuy(ipiv,jpiv) +  &
                 depthv(ipiv,jpiv)*flagv*facv*vbavg(ipiv,jpiv)*scvx(ipiv,jpiv) 
            ivoltrp(isec) = ivoltrp(isec)                                     +  &
                 volu(ipiv,jpiv)*flagu*facu*uice(ipiv,jpiv)*scuy(ipiv,jpiv) +  &
                 volv(ipiv,jpiv)*flagv*facv*vice(ipiv,jpiv)*scvx(ipiv,jpiv) 

            ! Negative transport - contributes if velocity X direction flag < 0
            facu=step(0.0,ubavg(ipiv,jpiv)*flagu)
            facv=step(0.0,vbavg(ipiv,jpiv)*flagv)
            voltrm(isec) = voltrm(isec)                                          +  &
                 depthu(ipiv,jpiv)*flagu*facu*ubavg(ipiv,jpiv)*scuy(ipiv,jpiv) +  &
                 depthv(ipiv,jpiv)*flagv*facv*vbavg(ipiv,jpiv)*scvx(ipiv,jpiv) 
            ivoltrm(isec) = ivoltrm(isec)                                     +  &
                 volu(ipiv,jpiv)*flagu*facu*uice(ipiv,jpiv)*scuy(ipiv,jpiv) +  &
                 volv(ipiv,jpiv)*flagv*facv*vice(ipiv,jpiv)*scvx(ipiv,jpiv) 
         end do ! points on section
    

         if (textfileflag) then 
            write(css,'(i3.3)') isec
            ! Write to file -- depth integrated transport
            fname='transport_net'//css//'.dat'
            if (.not.vapp) then
               open(33,file=trim(fname),status="replace",form="formatted",position="rewind")
               write(33,'(4a20)') '%  Time[year]','Vol trans[Sv]','Pos Vol Trans[Sv]', &
               'Neg Vol Trans[Sv]' 
            else
               open(33,file=trim(fname),status="unknown",form="formatted",position="append")
            end if
            write(33,'(f20.5,3e20.7)') rtime,voltr(isec)*1e-6,voltrp(isec)*1e-6, voltrm(isec)*1e-6
            close(33)
            write(345,'(a)') trim(fname)

            ! Write to file -- ice transport
            fname='transport_net_ice'//css//'.dat'
            if (.not.vapp) then
               open(33,file=trim(fname),status="replace",form="formatted",position="rewind")
               write(33,'(4a20)') '%  Time[year]','Vol trans[Sv]', 'Pos Vol Trans[Sv]', &
               'Neg Vol Trans[Sv]'
            else
               open(33,file=trim(fname),status="unknown",form="formatted",position="append")
            end if
            write(33,'(f20.5,3e20.7)') rtime,ivoltr(isec)*1e-6,ivoltrp(isec)*1e-6,ivoltrm(isec)*1e-6
            close(33)
            write(345,'(a)') trim(fname)

         end if


      end do ! section loop

      if (textfileflag) close(345)

      ! Special time case for calendar case
      if (hfile%yrflag == 3 .and. trim(hfile%ftype)/='nersc_weekly') then
         rtime2=nint(jday1950(hfile)*24)
         timeunits='hours since 1950-01-01 00:00:00'
      else
         rtime2=rtime
         timeunits=''
      end if

      ! KAL - new.. Put transport into netcdf file 
      call ncwrite_transportdata('transport'       ,voltr  *1e-6,nsec,rtime2,1, &
         appendfile=vapp,units='1e6 m3s-1',timeunits=trim(timeunits))
      call ncwrite_transportdata('transport_neg'   ,voltrm *1e-6,nsec,rtime2,1,&
         appendfile=vapp,units='1e6 m3s-1',timeunits=trim(timeunits))
      call ncwrite_transportdata('transport_pos'   ,voltrp *1e-6,nsec,rtime2,1,&
         appendfile=vapp,units='1e6 m3s-1',timeunits=trim(timeunits))
      call ncwrite_transportdata('icetransport'    ,ivoltr *1e-6,nsec,rtime2,1,&
         appendfile=vapp,units='1e6 m3s-1',timeunits=trim(timeunits))
      call ncwrite_transportdata('icetransport_neg',ivoltrm*1e-6,nsec,rtime2,1,&
         appendfile=vapp,units='1e6 m3s-1',timeunits=trim(timeunits))
      call ncwrite_transportdata('icetransport_pos',ivoltrp*1e-6,nsec,rtime2,1,&
         appendfile=vapp,units='1e6 m3s-1')
   end subroutine transport



   subroutine transport2(hfile,appendfile)
      use mod_grid
      use mod_sections
      use mod_hycomfile_io
      implicit none
      type(hycomfile), intent(in) :: hfile
      logical, intent(in), optional :: appendfile

      integer :: n,j,k,sec,ind,i,layer,ios,im1,jm1,ip1,jp1
      real, parameter :: cpsw = 3987.
      real, parameter :: t0deg    =273.15   ! 0 deg C in K
      real, parameter ::     g    =9.806    ! 0 deg C in K

      real, dimension(idm,jdm) :: olddepth,sumdepth, saln, dp, utot,vtot,temp, &
         transfacu,transfacv, salu, salv, temu, temv, dpu,dpv, pu, pv, &
         dens
      real, dimension(max_sdm) :: rngmsku, rngmskv, flwmsku, flwmskv
      character(len=3) :: css
      character(len=80) :: fname
      character(len=40) :: timeunits
      real, dimension(maxntrans) :: voltrans,heattrans ! NB - heattrans does not consider ice
      real, dimension(maxntrans) :: ivoltrans,iareatrans
      integer :: itrans, ipnt, isec, ipiv, jpiv, kdm ,icrit
      real :: masku, maskv,  rtime, rtime2
      logical :: vapp
      real :: thref
      include 'stmt_funcs0.H'

      ! set to true to use old text files
      logical, parameter :: textfileflag=.false. 

#if defined (SCALAR_TRANS) 
      integer :: ist
      character(len=80) :: cline, cline2, tmpstr
      character(len=2) :: cnst

      ! this sets up arrays to hold scalar transports, and name of fields to
      ! retrieve from datafile
      call scalar_transport_specification(idm,jdm)
#endif

      thref=1e-3

      vapp=.false.
      if (present(appendfile)) vapp=appendfile

      ! Start reading/processing layer by layer
      sumdepth=0.
      olddepth=0.
      voltrans(:)=0.
      heattrans(:)=0.
      pu=0.
      pv=0.
      kdm=vdim(hfile)
      rtime=fyear(hfile)
      do k=1,kdm

         !call HFReadDPField(hfile,dp,idm,jdm,k,1)
         call HFReadDPField_m(hfile,dp,idm,jdm,k,1) ! Always returns thickness in meters
         call HFReaduvtot  (hfile,utot,vtot,idm,jdm,k,1)
         call HFReadField  (hfile,temp,idm,jdm,'temp    ',k,1)
         if (trim(hfile%ftype)/='archv') then
            call HFReadField  (hfile,saln,idm,jdm,'saln    ',k,1) ! need "salin" for archv files
         else
            call HFReadField  (hfile,saln,idm,jdm,'salin   ',k,1) 
         end if

#if defined (SCALAR_TRANS) 
         ! Test code for calc. transport of a specified scalar variable (to be set in infile?)
         do i=1,num_scalar_trans
            call HFReadField  (hfile,scalar(:,:,i),idm,jdm,scalar_trans_name(i),k,1)
         end do
#endif

         ! Cumulative depth
         do j=1,jdm
         do i=1,idm
            olddepth(i,j)=sumdepth(i,j)
            sumdepth(i,j)= sumdepth(i,j) + dp(i,j)
         end do
         end do

         ! dp, sal, temp  in u/v points
         ! sal/temp/dens in u/v  ? Use upstream in stead?
         do j=1,jdm
         do i=1,idm
            if (periodic) then
               im1=mod(idm+i-2,idm)+1
            else
               im1=max(i-1,1)
            end if
            jm1=max(1,j-1)
            dpu(i,j)=max(0., &
                 min(depthu(i,j),.5*(sumdepth(i,j)+sumdepth(im1,j)))-&
                 min(depthu(i,j),.5*(olddepth(i,j)+olddepth(im1,j))))
            dpv(i,j)=max(0.,&
                 min(depthv(i,j),.5*(sumdepth(i,j)+sumdepth(i,jm1)))-&
                 min(depthv(i,j),.5*(olddepth(i,j)+olddepth(i,jm1))))
            pu(i,j)=min(depthu(i,j),.5*(olddepth(i,j)+olddepth(im1,j)))
            pv(i,j)=min(depthv(i,j),.5*(olddepth(i,j)+olddepth(i,jm1)))

            ! Salinities in u/v points
            salu(i,j)=dp(i,j)*saln(i,j)+dp(im1,j)*saln(im1,j)
            salv(i,j)=dp(i,j)*saln(i,j)+dp(i,jm1)*saln(i,jm1)
            salu(i,j) = 0.5 * salu(i,j) / (max(dpu(i,j),1e-4))
            salv(i,j) = 0.5 * salv(i,j) / (max(dpv(i,j),1e-4))
         
            ! Temperatures
            temu(i,j)=dp(i,j)*temp(i,j) + dp(im1,j)*temp(im1,j)
            temv(i,j)=dp(i,j)*temp(i,j) + dp(i,jm1)*temp(i,jm1)
            temu(i,j) = 0.5 * temu(i,j) /max(dpu(i,j),1e-4)
            temv(i,j) = 0.5 * temv(i,j) /max(dpv(i,j),1e-4)

            ! Densities
            dens(i,j)=sig(temp(i,j),saln(i,j))

         end do
         end do


         !  This is the total transport in this layer and the given grid cell
         do j=1,jdm
         do i=1,idm
            transfacu(i,j) = dpu(i,j)*utot(i,j)*scuy(i,j)
            transfacv(i,j) = dpv(i,j)*vtot(i,j)*scvx(i,j)
         end do
         end do


         ! Go through each section now
         do itrans=1,ntrans

            isec = transport_secnum(itrans)
            rngmsku(:)=1.
            rngmskv(:)=1.
            do icrit=1,transport_nbcrit(itrans)
               ! Set mask depending on variable range
               call rangemask(temu,temv,salu,salv,dpu,dpv,pu,pv,  &
                  itrans,rngmsku,rngmskv,isec,icrit)
            enddo
            ! Set mask depending on flow direction
            call flowmask(utot, vtot, flwmsku, flwmskv,isec, itrans)

            ! Set final mask - combines flow direction mask and range mask
            do ipnt=1,sdm(isec)
               ipiv=ndeipiv(ipnt,isec)
               jpiv=ndejpiv(ipnt,isec)
               
               masku = flwmsku(ipnt) * rngmsku(ipnt) * ndeflagu(ipnt,isec)
               maskv = flwmskv(ipnt) * rngmskv(ipnt) * ndeflagv(ipnt,isec)

               voltrans(itrans) = voltrans(itrans) +  &
                    masku*transfacu(ipiv,jpiv)     +  &
                    maskv*transfacv(ipiv,jpiv) 

               !  NB -- requires dp to be in pressure coords - is relative to 0
               !  deg C
               heattrans(itrans) = heattrans(itrans)                                           +  &
                    masku*transfacu(ipiv,jpiv)*((temp(ipiv,jpiv))-0.)*cpsw*(1000.) + &  !+dens(ipiv,jpiv))  +  &
                    maskv*transfacv(ipiv,jpiv)*((temp(ipiv,jpiv))-0.)*cpsw*(1000.)  !+dens(ipiv,jpiv))
#if defined (SCALAR_TRANS) 
               do ist=1,num_scalar_trans
                  scalar_trans(itrans,ist) = scalar_trans(itrans,ist) +   &
                    scalar_trans_factor(ist)* ( &
                       masku*transfacu(ipiv,jpiv)*(scalar(ipiv,jpiv,ist) - scalar_trans_offset(ist)) + &
                       maskv*transfacv(ipiv,jpiv)*(scalar(ipiv,jpiv,ist) - scalar_trans_offset(ist))   &
                    )
               end do
#endif
            end do
         end do !itrans
      end do !kdm


      if (textfileflag) then 
      open(345,file='section_transport2.filelist',status='replace')
      do itrans=1,ntrans

         voltrans(itrans)=voltrans(itrans)
         sec = transport_secnum(itrans)

         print '(a,i3,a,a21,a,i1)', &
            'Transport ',itrans,' With ID ',             &
            transport_name(itrans),'# criterium',transport_nbcrit(itrans)
            do icrit=1,transport_nbcrit(itrans)
               print '(a,a14,2f10.2)',' range ',transport_rangeid(itrans,icrit), &
                  transport_lowrange(itrans,icrit), transport_uprange(itrans,icrit)
            enddo
            
            
          print '(a,a1,a,i3,a,f16.4)', &
            ' Direction ',transport_directionid(itrans), &
            ' section ',sec,  &
            ' transport:',voltrans(itrans) * 1e-6 


         write(css,'(i3.3)') sec
 
         ! Write to file -- depth integrated transport
         fname=trim(transport_name(itrans))//'.dat'
         if (.not.vapp) then
            open(33,file=trim(fname),status="replace",form="formatted",position="rewind")
            write (33,'(a)') trim(transport_name(itrans)) 
            do icrit=1,transport_nbcrit(itrans)
               write (33,'(a,a14,2f10.2)') &
               ' with range type and limits ',transport_rangeid(itrans,icrit), &
               transport_lowrange(itrans,icrit), transport_uprange(itrans,icrit)           
            enddo
            write (33,'(a,a1,a,a20)') &
               ' with direction ',transport_directionid(itrans), &
               ' across section ',trim(secname(sec))
#if defined (SCALAR_TRANS) 
            tmpstr=''
            if (num_scalar_trans>0) then
               tmpstr='  +  Scalar transport of: '
               do ist=1,num_scalar_trans
                  tmpstr=trim(tmpstr)//trim(scalar_trans_name(ist)(:))//','
               end do
               tmpstr=tmpstr(1:len_trim(tmpstr)-1)
               write (33,'(a)') 'Columns are time [year] and transport [Sv], Heat transport[W]'//trim(tmpstr)
            else
               write (33,'(a)') 'Columns are time [year] and transport [Sv], Heat transport[W]'
            end if
#else
            write (33,'(a)') 'Columns are time [year] and transport [Sv], Heat transport[W]'
#endif
         else
            open(33,file=trim(fname),status="unknown",form="formatted",position="append")
         end if
         !open(33,file=trim(fname),status="replace",form="formatted",position="rewind")
#if defined (SCALAR_TRANS) 
         write(cnst,'(i2)') num_scalar_trans
         !write(cline,'(f14.4,f14.4,e14.4)') rtime,voltrans(itrans)*1e-6,heattrans(itrans) !MS. 120312
         write(cline,'(f14.4,f18.4,e14.4)') rtime,voltrans(itrans)*1e-6,heattrans(itrans)
         cline2=''
         if (num_scalar_trans>0) then
            write(cline2,'('//cnst//'e14.4)') scalar_trans(itrans,:)
         end if
         write(33,'(a)') trim(cline)//trim(cline2)
#else
         write(33,'(f14.4,f14.4,e14.4)') rtime,voltrans(itrans)*1e-6,heattrans(itrans)
#endif
         close(33)


         ! This is a file list generated for each pass - dumps the transport
         ! files which are created

         write(345,'(a)') trim(transport_name(itrans))

      end do
      close(345)
      end if !(textfileflag) 

      ! Special time case for calendar case
      if (hfile%yrflag == 3 .and. trim(hfile%ftype)/='nersc_weekly') then
         rtime2=nint(jday1950(hfile)*24)
         timeunits='hours since 1950-01-01 00:00:00'
      else
         rtime2=rtime
         timeunits=''
      end if


      call ncwrite_transportdata('vol_transport',voltrans(1:ntrans)*1e-6,ntrans,rtime2,2,units='m3 s-1', &
         comment='Volume Transport',appendfile=vapp,timeunits=trim(timeunits))
      call ncwrite_transportdata('heat_transport',heattrans(1:ntrans),ntrans,rtime2,2,units='W', &
         comment='Heat Transport',appendfile=vapp,timeunits=trim(timeunits))
#if defined (SCALAR_TRANS) 
      do ist=1,num_scalar_trans
         call ncwrite_transportdata(trim(scalar_trans_name(ist))//'_transport', &
            scalar_trans(:,ist),ntrans,rtime2,2,appendfile=vapp,timeunits=trim(timeunits))
      end do
#endif
   end subroutine transport2



   subroutine ncwrite_transportdata(cfld,trnsp,nsecin,rtime,mode,appendfile,&
      units,comment,timeunits)
   use m_handle_err
   use netcdf
   use mod_sections
   implicit none
   integer,          intent(in) :: nsecin, mode
   character(len=*), intent(in) :: cfld
   real,             intent(in) :: trnsp(nsecin), rtime
   logical, intent(in), optional :: appendfile
   character(len=*), intent(in), optional :: units, comment, timeunits

   character(len=80) :: ncfil
   logical :: ex,v3d, vtime, vapp
   integer :: ncid, varid, ierr, isec,v2d(2), &
      sect_dimension_id, time_dimension_id, c20id, itrans, &
      varidtrname, varidsecname, varidmask, lonid, latid, n2id, &
      varidmaskrange, dirid, timevarid
   integer, save :: rdimlen=-1


   vapp=.false.
   if (present(appendfile)) vapp=appendfile


   ! Mode is either 1 or 2. Mode 1 is used by transport, mode 2 by transport2.
   ! The main difference is that mode 1 loops over sections, whereas mode 2
   ! loops over defined transports
   if (mode==1) then
      if (nsecin/=nsec) then
         print *,'ncwrite_transportdata - mode 1:'
         print *,'input nsec must match mod_sections nsec'
         stop '(mod_transport:ncwrite_transportdata)'
      end if
      ncfil = 'transports.nc'
   else if (mode==2) then
      if (nsecin/=ntrans) then
         print *,'ncwrite_transportdata - mode 2:'
         print *,'input nsec must match mod_transport ntrans'
         stop '(mod_transport:ncwrite_transportdata)'
      end if
      ncfil = 'transports2.nc'
   else
      print *,'ncwrite_transportdata - mode:'
      print *,'input mode must be 1 or 2 '
      stop '(mod_transport:ncwrite_transportdata)'
   end if
   inquire(exist=ex,file=trim(ncfil))

   ! Base definition of file if it does not exist, or 
   ! if appendfile is switched off AND this is the first call
   if (.not. ex .or. (.not. vapp .and. rdimlen==-1)) then
      !print *,'NC DEFINE!'
      ! Define dimensions
      call handle_err(NF90_create(trim(ncfil),NF90_CLOBBER,ncid))
      call handle_err(NF90_DEF_DIM(ncid,'time',NF90_UNLIMITED,time_dimension_id))
      call handle_err(NF90_DEF_VAR(ncid,'time',NF90_FLOAT,time_dimension_id,timevarid))
      if (present(timeunits)) then
      if (trim(timeunits)/='') then
         call handle_err(NF90_PUT_ATT(ncid,timevarid , 'units',timeunits))
      end if
      end if
      call handle_err(NF90_DEF_DIM(ncid,'char20' ,20,c20id))
      call handle_err(NF90_DEF_DIM(ncid,'N_2',2,N2id))

      ! Define global attributes
      call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'title','TOPAZ transport across sections'))
      call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'institution', &
         'NERSC, Thormoehlens gate 47, N-5006 Bergen, Norway'))
      call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'references', &
         'http://topaz.nersc.no'))
      call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'comment', &
      'Transports are defined as positive if it is towards the right-hand side &
      when standing at the start of a section and looking towards the end of the &
      section'))


      rdimlen=1 
      if (mode==1) then ! barotropic transport
         call handle_err(NF90_DEF_DIM(ncid,'section',nsec,sect_dimension_id))
         call handle_err(NF90_DEF_VAR(ncid,'sectionname' , NF90_CHAR, (/c20id,sect_dimension_id/),varid))
         call handle_err(NF90_PUT_ATT(ncid,varid , 'comment','Name of section'))
         call handle_err(NF90_DEF_VAR(ncid,'section_lon' , NF90_FLOAT, (/N2id,sect_dimension_id/),lonid))
         call handle_err(NF90_PUT_ATT(ncid,lonid , 'comment','Start and end longitudes of section '))
         call handle_err(NF90_DEF_VAR(ncid,'section_lat' , NF90_FLOAT, (/N2id,sect_dimension_id/),latid))
         call handle_err(NF90_PUT_ATT(ncid,latid , 'comment','Start and end latitudes of section '))
         call handle_err(nf90_enddef(ncid))
         do isec=1,nsec 
            call handle_err(nf90_put_var(ncid,varid,secname(isec)(1:20),(/1,isec/)))
            call handle_err(nf90_put_var(ncid,lonid ,(/ndelon(1,isec),ndelon(sdm(isec),isec)/),(/1,isec/)))
            call handle_err(nf90_put_var(ncid,latid ,(/ndelat(1,isec),ndelat(sdm(isec),isec)/),(/1,isec/)))
         end do
      else if (mode==2) then ! detailed transport
         call handle_err(NF90_DEF_DIM(ncid,'transport',ntrans,sect_dimension_id))
         call handle_err(NF90_DEF_VAR(ncid,'sectionname' , NF90_CHAR, (/c20id,sect_dimension_id/),varidsecname))
         call handle_err(NF90_PUT_ATT(ncid,varidsecname , 'comment','Name of section'))
         call handle_err(NF90_DEF_VAR(ncid,'transportname' , NF90_CHAR, (/c20id,sect_dimension_id/),varidtrname))
         call handle_err(NF90_DEF_VAR(ncid,'mask_variable_name' , NF90_CHAR, (/c20id,sect_dimension_id/),varidmask))
         call handle_err(NF90_PUT_ATT(ncid,varidmask , 'comment','Mask variable for transport'))
         call handle_err(NF90_DEF_VAR(ncid,'mask_variable_range' , NF90_FLOAT, (/N2id,sect_dimension_id/),varidmaskrange))
         call handle_err(NF90_PUT_ATT(ncid,varidmaskrange , 'comment','Mask variable range for transport'))
         call handle_err(NF90_DEF_VAR(ncid,'section_lon' , NF90_FLOAT, (/N2id,sect_dimension_id/),lonid))
         call handle_err(NF90_PUT_ATT(ncid,lonid , 'comment','Start and end longitudes of section '))
         call handle_err(NF90_DEF_VAR(ncid,'section_lat' , NF90_FLOAT, (/N2id,sect_dimension_id/),latid))
         call handle_err(NF90_PUT_ATT(ncid,latid , 'comment','Start and end latitudes of section '))
         call handle_err(NF90_DEF_VAR(ncid,'direction' , NF90_CHAR, (/sect_dimension_id/),dirid))
         call handle_err(NF90_PUT_ATT(ncid,dirid , 'comment','direction flag + is positive, - is negative, n is net'))
         call handle_err(nf90_enddef(ncid))
         !print *,'hei3'
         do itrans=1,ntrans 
            isec=transport_secnum(itrans)
            call handle_err(nf90_put_var(ncid,varidsecname,secname(isec)(1:20),(/1,itrans/)))
            call handle_err(nf90_put_var(ncid,varidtrname ,transport_name(itrans)(1:20),(/1,itrans/)))
            if (transport_nbcrit(itrans).gt.1) then
               call handle_err(nf90_put_var(ncid,varidmask ,'MULTIPLE',(/1,itrans/)))
               call handle_err(nf90_put_var(ncid,varidmaskrange ,(/9999, &
                                         9999/),(/1,itrans/)))
            else
               call handle_err(nf90_put_var(ncid,varidmask ,transport_rangeid(itrans,1)(1:12),(/1,itrans/)))
               call handle_err(nf90_put_var(ncid,varidmaskrange ,(/transport_lowrange(itrans,1), &
                                         transport_uprange(itrans,1)/),(/1,itrans/)))
            end if
            call handle_err(nf90_put_var(ncid,lonid ,(/ndelon(1,isec),ndelon(sdm(isec),isec)/),(/1,itrans/)))
            call handle_err(nf90_put_var(ncid,latid ,(/ndelat(1,isec),ndelat(sdm(isec),isec)/),(/1,itrans/)))
            call handle_err(nf90_put_var(ncid,dirid ,transport_directionid(itrans), (/itrans/)))
         end do
      end if
   else
      !print *,'NC RETRIEVE!'
      ! Get dimensions - we may need them for defining variables
      call handle_err(NF90_open(trim(ncfil),NF90_WRITE,ncid))
      call handle_err(nf90_inq_dimid(ncid, 'time'  , time_dimension_id))
      if (mode==1) then ! barotropic transport
         call handle_err(NF90_inq_dimid(ncid,'section'  ,sect_dimension_id))
      else if (mode==2) then ! detailed transport
         call handle_err(NF90_inq_dimid(ncid,'transport'  ,sect_dimension_id))
      end if
      if (rdimlen == -1 ) then
         call handle_err(NF90_INQUIRE_DIMENSION(ncid,time_dimension_id,len=rdimlen))
         rdimlen=rdimlen+1
      end if
      call handle_err(nf90_inq_varid(ncid,'time',timevarid))
   end if
   v2d=(/sect_dimension_id,time_dimension_id/)



   !print *,'NC WRITE!'
   call handle_err(nf90_put_var(ncid,timevarid,rtime,start=(/rdimlen/)))
   if (mode==1) then ! barotropic transport

      ! Inquire and define (if necessary) and put the variable
      ierr=nf90_inq_varid(ncid,trim(cfld),varid)
      if (ierr/=NF90_NOERR) then
         call handle_err(nf90_redef(ncid))
         call handle_err(nf90_def_var(ncid,trim(cfld),NF90_Float,v2d,varid))
         if (present(units)) call handle_err(nf90_put_att(ncid, varid, 'unit', units))
         if (present(comment)) call handle_err(nf90_put_att(ncid, varid, 'comment', units))
         call handle_err(nf90_enddef(ncid))
      end if
      call handle_err(nf90_put_var(ncid,varid,trnsp(1:nsec),start=(/1,rdimlen/)))
      call handle_err(nf90_close(ncid))


   else  ! detailed transports

      ! Inquire and define (if necessary) and put the variable
      ierr=nf90_inq_varid(ncid,trim(cfld),varid)
      if (ierr/=NF90_NOERR) then
         call handle_err(nf90_redef(ncid))
         call handle_err(nf90_def_var(ncid,trim(cfld), NF90_Float,v2d,varid))
         if (present(units)) call handle_err(nf90_put_att(ncid, varid, 'unit', units))
         if (present(comment)) call handle_err(nf90_put_att(ncid, varid, 'comment', units))
         call handle_err(nf90_enddef(ncid))
      end if
      call handle_err(nf90_put_var(ncid,varid,trnsp(1:ntrans),start=(/1,rdimlen/)))
      call handle_err(nf90_close(ncid))
   end if
   end subroutine



   real function step(xx1,xx2)
      real,intent(in) :: xx1,xx2
      step=(1.+sign(1.,xx1-xx2))*.5
   end function step


   subroutine rangemask(temu,temv,salu,salv,dpu,dpv,pu,pv, &
      itrans,secmsku,secmskv,isec,icrit)
   use mod_xc
   use mod_sections
   use m_depth_frac
   implicit none
   real, dimension(idm,jdm), intent(in) :: temu, temv, salu, salv,  &
      dpu,dpv,pu,pv
   integer, intent(in) :: itrans, isec, icrit
   real,   intent(inout), dimension(max_sdm) :: secmsku,secmskv
   integer :: ipiv, jpiv, ipnt
   real,   dimension(max_sdm) :: maskvaru, maskvarv, masku, maskv
   real :: dfracu(idm,jdm),dfracv(idm,jdm)
   real :: pi,radian, thref
   include 'stmt_funcs0.H'
   thref=1e-3
   pi=atan(1.)*4
   radian=180./pi

   ! Range type
   if (transport_rangeid(itrans,icrit)=='SALINITY') then
      do ipnt=1,sdm(isec)
         ipiv=ndeipiv(ipnt,isec)
         jpiv=ndejpiv(ipnt,isec)
         maskvaru(ipnt)=salu(ipiv,jpiv)
         maskvarv(ipnt)=salv(ipiv,jpiv)
      end do
   else if (transport_rangeid(itrans,icrit)=='TEMPERATURE') then
      do ipnt=1,sdm(isec)
         ipiv=ndeipiv(ipnt,isec)
         jpiv=ndejpiv(ipnt,isec)
         maskvaru(ipnt)=temu(ipiv,jpiv)
         maskvarv(ipnt)=temv(ipiv,jpiv)
      end do
   else if (transport_rangeid(itrans,icrit)=='DENSITY') then
      do ipnt=1,sdm(isec)
         ipiv=ndeipiv(ipnt,isec)
         jpiv=ndejpiv(ipnt,isec)
         maskvaru(ipnt)=sig(temu(ipiv,jpiv),salu(ipiv,jpiv))
         maskvarv(ipnt)=sig(temv(ipiv,jpiv),salv(ipiv,jpiv))
      end do
   else if (transport_rangeid(itrans,icrit)=='NONE') then
      do ipnt=1,sdm(isec)
         masku(ipnt)=1.
         maskv(ipnt)=1.
      end do
   else if (transport_rangeid(itrans,icrit)/='DEPTH') then
      print *,'Unknown range type '//&
      transport_rangeid(itrans,icrit)
      call exit(1)
   end if
   ! Set range mask
   if (transport_rangeid(itrans,icrit)=='SALINITY'    .or. &
       transport_rangeid(itrans,icrit)=='TEMPERATURE' .or. &
       transport_rangeid(itrans,icrit)=='DENSITY') then
      do ipnt=1,sdm(isec)
         ipiv=ndeipiv(ipnt,isec)
         jpiv=ndejpiv(ipnt,isec)
         masku(ipnt) = secmsku(ipnt)*step(maskvaru(ipnt),transport_lowrange(itrans,icrit))
         masku(ipnt) = masku(ipnt) * step(transport_uprange(itrans,icrit),maskvaru(ipnt))
         maskv(ipnt)= secmskv(ipnt)*step(maskvarv(ipnt),transport_lowrange(itrans,icrit))
         maskv(ipnt)= maskv(ipnt) * step(transport_uprange(itrans,icrit),maskvarv(ipnt))
      end do
   elseif  (transport_rangeid(itrans,icrit)=='DEPTH' ) then
      ! NB - lowrange/uprange refers to values, limits in depth_frac
      ! refer to depth, with positive dir downward. So uprange is
      ! actually depth of lowest layer
      call depth_frac(transport_lowrange(itrans,icrit),transport_uprange(itrans,icrit),pu,pu+dpu,dfracu,idm,jdm)
      call depth_frac(transport_lowrange(itrans,icrit),transport_uprange(itrans,icrit),pv,pv+dpv,dfracv,idm,jdm)
      do ipnt=1,sdm(isec)
         ipiv=ndeipiv(ipnt,isec)
         jpiv=ndejpiv(ipnt,isec)
         masku(ipnt)=secmsku(ipnt)*dfracu(ipiv,jpiv)
         maskv(ipnt)=secmskv(ipnt)*dfracv(ipiv,jpiv)
      end do
   end if

   do ipnt=1,sdm(isec)
      secmsku(ipnt)=masku(ipnt)
      secmskv(ipnt)=maskv(ipnt)
   end do
   end subroutine


   subroutine flowmask(utot, vtot, flwmsku, flwmskv,isec, itrans)
   use mod_xc
   use mod_sections
   implicit none
   integer, intent(in)  :: isec, itrans
   real, dimension(idm,jdm), intent(in)   :: utot, vtot
   real, dimension(max_sdm), intent(out)  :: flwmsku, flwmskv
   real  :: steplus, stepminus, stepnet, stepv, stepu, &
      omstepu, omstepv, flagu, flagv
   integer :: ipiv, jpiv, ipnt

   ! Set direction mask
   steplus=0. ; stepminus=0.; stepnet = 0.
   if(transport_directionid(itrans)=='+') steplus   =1.
   if(transport_directionid(itrans)=='-') stepminus =1.
   if(transport_directionid(itrans)=='n') stepnet   =1.

   if(abs(steplus+stepminus+stepnet)<1e-4) then
      print *,'unknown direction id for transport '// &
         trim(transport_name(itrans))
      call exit(1)
   end if

   ! Set step factors
   do ipnt=1,sdm(isec)
      
      ipiv=ndeipiv(ipnt,isec)
      jpiv=ndejpiv(ipnt,isec)
      flagu=ndeflagu(ipnt,isec)
      flagv=ndeflagv(ipnt,isec)

      stepu  =step( utot(ipiv,jpiv)*flagu,0.) ! U Transp > 0. => 1
      stepv  =step( vtot(ipiv,jpiv)*flagv,0.) ! V Transp > 0. => 1
      omstepu=step(0.0,utot(ipiv,jpiv)*flagu)       ! U Transp < 0. => 1
      omstepv=step(0.0,vtot(ipiv,jpiv)*flagv)       ! V transp < 0. => 1

      flwmsku(ipnt)=stepu * steplus ! = 1 if u transp > 0 & dir = +
      flwmskv(ipnt)=stepv * steplus ! = 1 if v transp > 0 & dir = +

      flwmsku(ipnt)=max(omstepu * stepminus,flwmsku(ipnt)) ! = 1 if u transp < 0 & dir = - or previous
      flwmskv(ipnt)=max(omstepv * stepminus,flwmskv(ipnt)) ! = 1 if v transp < 0 & dir = - or previous

      flwmsku(ipnt)=max(stepnet,flwmsku(ipnt)) ! = 1 if net u transp, or previous
      flwmskv(ipnt)=max(stepnet,flwmskv(ipnt)) ! = 1 if net v transp, or previous

   end do
   end subroutine

end module mod_transport

