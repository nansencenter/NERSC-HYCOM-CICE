      program clim_weekly
      use mod_xc
      use mod_za
      use mod_grid
      use mod_year_info
      use mod_hycomfile_io
      use m_parse_blkdat
      implicit none

      character(len= *),parameter ::fnestin ='clim_weekly.in'

      character(len=80) :: cfile,cname,c80, ftype
      character(len=80) :: filelist(200)
      character(len=10) :: filetype
      character(len=3)  :: rungen,oldrungen
      character(len=2)  :: cmonth
      character(len=1)  :: cweek

      ! Arrays holding data on inner grid (where we need nesting
      ! conditions)
      real, allocatable, dimension(:,:) ::  &
         temp,saln,dp, tmpdp, tmpfld, &
         hicem,ficem,ssh,hsnwm,uice,vice,utot,vtot
      type(year_info)   :: rt
      integer, allocatable, dimension(:,:) ::  &
         idummy

      real :: realvar,amin,amax
      integer :: i,j,k,ios,klevel, irec_offset,ivar, &
         kdm, intvar,indx
         integer :: imonth, iweek, imonth2, iweek2, itime, nfiles,ifile
      logical :: ex, isvelocity,ass

      ! Only set to true for testing output with "plotfp". If set to true,
      ! the routines working with weekly average files WILL NOT WORK
      logical :: testoutput=.false. 

      type(hycomfile) :: hfile



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set up global  grids ...  There is redundant information here, blkdat.input
! contains idm,jdm,kdm. regional.depths.b contains idm jdm

      ! 1 - start with xcspmd, which accesses regional.grid.b
      ! Initialize Arrai IO
      call xcspmd()
      call zaiost()
      call get_grid()

      ! 2 - get kdm from blkdat.input - cross check idm/kdm with existing
      ! Search for idm/jdm/kdm
      print *,'parsing blkdat (consistency check  + kdm)'
      
      call parse_blkdat('idm   ','integer',realvar,intvar)
      if (intvar/=idm) then
         print *,'idm: blkdat.input and regional.depths.b mismatch'
         stop '(nest_offline)'
      end if
      call parse_blkdat('jdm   ','integer',realvar,intvar)
      if (intvar/=jdm) then
         print *,'jdm: blkdat.input and regional.depths.b mismatch'
         stop '(nest_offline)'
      end if
      call parse_blkdat('kdm   ','integer',realvar,intvar)
      kdm=intvar




      ! Allocate variables for global grid
      allocate(idummy(idm,jdm)) ; idummy=0
      allocate(temp  (idm,jdm))
      allocate(saln  (idm,jdm))
      allocate(dp    (idm,jdm))
      allocate(utot  (idm,jdm))
      allocate(vtot  (idm,jdm))
      allocate(tmpfld(idm,jdm))
      allocate(tmpdp (idm,jdm))

      allocate(uice(idm,jdm))
      allocate(vice(idm,jdm))
      allocate(hicem(idm,jdm))
      allocate(hsnwm(idm,jdm))
      allocate(ficem(idm,jdm))
      allocate(  ssh(idm,jdm))





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! We are set .. Now we must read the weekly average files from the
! global model, calculate monthly/weekly statistics, then dump new
! weekly files.

      inquire(exist=ex,file=fnestin)
      ! Contains weekly average files, one per line. They must be sorted
      ! so that all week1 month1 files are present first, all week2,
      ! month2 files are second and so on.


      if (.not. ex) then
         print *,'File "'//fnestin//' is not present'
         stop '(nest_climatology_from_weekly)'
      else
         open(118,file=fnestin)
      end if



      ! Get files to read from nestoffline.in
      read(118,*,iostat=ios) cfile
      oldrungen=''
      do itime=1,48 ! 12 weeks, 48 days

         ! A comment sign (#) separate file collections
         if (cfile(1:1) /= "#") then
            print *,'input file error '
            stop
         end if

         ! The entries are supposed to be for this week, this month
         iweek =mod(itime-1,4)+1
         imonth=ceiling(itime/4.);

         ! Start reading entries
         ios=0
         nfiles=0
         do while( ios==0)
         
            read(118,*,iostat=ios) cfile


            if (cfile(1:1) /= "#") then ! assume valid entry

               ! From file name - retrieve month/week info
               read(cfile(13:14),*) imonth2
               read(cfile(16:16),*) iweek2

               if (iweek/=iweek2 .or. imonth/=imonth2) then
                  print *,'week/month mismatch'
                  print *,iweek,iweek2
                  print *,imonth,imonth2
                  stop
               end if

               !! Check filetype - only ave is accepted
               if (trim(cfile(4:6))=='AVE') then
                  filetype='weekly'
               else
                  print *,'Unknown filetype '//trim(cfile)
                  stop '(nest_offline)'
               end if


               ! Augment nfiles counter
               nfiles=nfiles+1

               ! Add to file list
               indx=max(index(cfile,'.a'),index(cfile,'b'))
               filelist(nfiles)=cfile(1:indx-1)

            else ! Assume next set of entries
               ios=1
            end if


         end do

         if (nfiles==0) then
            print *,'No files found for ',iweek, imonth
            stop
         else
            print *,nfiles,' files found for ',iweek, imonth
         end if


         rungen=filelist(1)(1:3)
         if (trim(oldrungen)=='') then
            oldrungen=rungen
         end if



         ! Open Climatology file
         write(cmonth,'(i2.2)') imonth
         write(cweek ,'(i1.1)') iweek 
         cname=rungen//'AVE_9999_'//cmonth//'_'//cweek
         call zaiopf(trim(cname)//'.a','replace',218)
         open(218,file=trim(cname)//'.b',status='replace',form='formatted')
         print '(a,i4,a,i4)','Processing for month ',imonth,' and week ',iweek

         ! Copy first in filelist ave header to the new file
         open(219,file=trim(filelist(1))//'.b',status='old')
         do i=1,12
            read (219,'(a)') c80
            write(218,'(a)') c80
         end do
         close(219)

         ! Init hycomfile info
         ftype=getfiletype(trim(cfile)) ! should be weekly
         call initHF(hfile,trim(cfile),trim(ftype))


         ! Next line is "ave" counter - set to nfiles now
         write(218,115) nfiles
         write(218,116) 


         ! Calc stuff goes here !
         do klevel=1,kdm

            ! 3D vars
            saln=0.
            temp=0.
            utot=0.
            vtot=0.
            dp  =0.

            ! 2D vars
            uice=0.
            vice=0.
            hicem=0.
            ficem=0.
            hsnwm=0.
            ssh=0.

            print *,'   - layer ',klevel
            do ifile=1,nfiles

               if (oldrungen/=rungen) then
                  print *,'Files have different rungen'
                  print *,'old rungen:',oldrungen
                  print *,'new rungen:',rungen
                  stop
               end if


               ! Process dp
               !call HFReadDPField(hfile,tmpdp,idm,jdm,klevel,1)
               call HFReadDPField_p(hfile,tmpdp,idm,jdm,klevel,1) ! Returns in pressure coords

               dp=dp+tmpdp

               ! Process salinities
               call HFReadField(hfile,tmpfld,idm,jdm,'saln    ',klevel,1)
               saln=saln+tmpfld

               ! Process salinities
               call HFReadField(hfile,tmpfld,idm,jdm,'temp    ',klevel,1)
               temp=temp+tmpfld

               ! Process utot      
               call HFReadField(hfile,tmpfld,idm,jdm,'utot    ',klevel,1)
               utot=utot+tmpfld

               ! Process vtot      
               call HFReadField(hfile,tmpfld,idm,jdm,'vtot    ',klevel,1)
               utot=utot+tmpfld


               ! ice variables
               if (klevel==1) then
                  ! Process ssh       
                  call HFReadField(hfile,tmpfld,idm,jdm,'ssh     ',0,1)
                  ssh=ssh+tmpfld

                  ! Process ficem      
                  call HFReadField(hfile,tmpfld,idm,jdm,'ficem   ',0,1)
                  ficem=ficem+tmpfld

                  ! Process hicem      
                  call HFReadField(hfile,tmpfld,idm,jdm,'hicem   ',0,1)
                  hicem=hicem+tmpfld

                  !! Process hsnwm      
                  !call read_weekfield2d(trim(filelist(ifile)),'hsnwm   ',tmpfld, &
                  !   idm,jdm,0,undef,lsilent=.true.)
                  !hsnwm=hsnwm+tmpfld

                  ! Process uice      
                  call HFReadField(hfile,tmpfld,idm,jdm,'uice    ',0,1)
                  uice =uice +tmpfld

                  ! Process vice      
                  call HFReadField(hfile,tmpfld,idm,jdm,'vice    ',0,1)
                  vice =vice +tmpfld
               end if

            end do

            if (klevel==1) then
               print *,'   - processing 2D variables '

               ! Makes it easier to verify variables have been read correctly
               if (testoutput) then
                  ssh=ssh/nfiles
                  ficem=ficem/nfiles
                  hicem=hicem/nfiles
                  hsnwm=hsnwm/nfiles
                  uice =uice /nfiles
                  vice =vice /nfiles
               end if

               call zaiowr(ssh  ,idummy,.false.,amin,amax,218,.false.)
               write(218,117) 'ssh     ',0,0.,0,0.,amin,amax 

               call zaiowr(ficem,idummy,.false.,amin,amax,218,.false.)
               write(218,117) 'ficem   ',0,0.,0,0.,amin,amax 

               call zaiowr(hicem,idummy,.false.,amin,amax,218,.false.)
               write(218,117) 'hicem   ',0,0.,0,0.,amin,amax 

               !call zaiowr(hsnwm,idummy,.false.,amin,amax,218,.false.)
               !write(218,117) 'hsnwm   ',0,0.,0,0.,amin,amax 

               call zaiowr(uice ,idummy,.false.,amin,amax,218,.false.)
               write(218,117) 'uice    ',0,0.,0,0.,amin,amax 

               call zaiowr(vice ,idummy,.false.,amin,amax,218,.false.)
               write(218,117) 'vice    ',0,0.,0,0.,amin,amax 
            end if

               if (testoutput) then
                  where (dp>1e-4)
                     utot=utot/dp
                     vtot=vtot/dp
                     saln=saln/dp
                     temp=temp/dp
                     dp  =dp/nfiles
                  elsewhere
                     utot=0.
                     vtot=0.
                     saln=0.
                     temp=0.
                     dp  =0.
                  end where
               end if

            ! The average values
            call zaiowr(dp  ,idummy,.false.,amin,amax,218,.false.)
            write(218,117) 'pres    ',0,0.,klevel,0.,amin,amax 

            call zaiowr(saln,idummy,.false.,amin,amax,218,.false.)
            write(218,117) 'saln    ',0,0.,klevel,0.,amin,amax 

            call zaiowr(temp,idummy,.false.,amin,amax,218,.false.)
            write(218,117) 'temp    ',0,0.,klevel,0.,amin,amax 

            call zaiowr(utot,idummy,.false.,amin,amax,218,.false.)
            write(218,117) 'utot    ',0,0.,klevel,0.,amin,amax 

            call zaiowr(vtot,idummy,.false.,amin,amax,218,.false.)
            write(218,117) 'vtot    ',0,0.,klevel,0.,amin,amax 

         end do
         call zaiocl(218)




      end do ! itime loop
 115  format (i5,4x,'''count '' = Averaging counter      ')
 116  format ( 'field       time step  model day', &
       '  k  dens        min              max')
 117  format (a8,' = ',i11,f11.2,i3,f7.3,1p2e16.7)

      end program clim_weekly
