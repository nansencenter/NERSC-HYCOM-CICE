      program multiproc
      use mod_xc
      use mod_za
      use mod_grid
      use mod_year_info
      use mod_pakk
      use m_pakk_field_ids
      implicit none

      type(year_info) rt

      character(len=40) outname
      character(len=3) rungen
      character(len=3) cmonth
      character(len=5) rforce
      character(len=5) pakk_id
      character(len=4) cyear
      character(len=11) tag7

      integer, parameter :: maxfiles=1000
      integer, parameter :: maxfields=1000


      integer nrfiles
      character(len=80) c80,dfile
      character(len=8 ) c8,actual_field
      character(len=80) fieldhead
      character(len=2 ), dimension(:), allocatable ::  util


      integer :: n2d,n3d
      character(len=5 ) char2d (maxfields)
      character(len=5 ) char3d (maxfields)
      character(len=5 ) varlist(maxfields)

      character(len=40) filenames (maxfiles)
      integer           fileunits (maxfiles)
      integer           avecounter(maxfiles)
      character(len=40) fileb
      integer piv,ip1,jp1
      real fac,rho
      integer ipiv,nrrec
      integer month,month1,month2,year
      logical ex,exq
      integer i,j,k,l,find,findab,finduf

      real,    dimension(:,:), allocatable ::
     &   h_acc,h, tmp, ut, vt, dp, tem, sal, speed,eke,mke,
     &   ub,vb,dpsum,ficem,hicem, fice, hice
      real :: sigma(100)

      real, parameter :: thref=1e-3
      real, parameter :: qthref=1e3
      real, parameter :: onem=9806.

      integer :: counter_acc

      integer :: nstep,klevel,actual_level
      real    :: coord
      real*8  :: time
      real*8, allocatable :: meanssh(:,:)

      integer :: ios, length, irec,maxlevel, idat(3)
      real    :: hmin, hmax, hmin2, hmax2
      integer :: sallevel,temlevel, dplevel, utlevel, vtlevel, splevel
      integer :: ficelevel, hicelevel
      logical :: goticestate

      !include 'stmt_fns.h'


      lp=6
      sigma=1.
      n2d=0
      n3d=0
      irec=1


      ! init IO
      call xcspmd()
      call zaiost()
      call get_grid()

      allocate(h_acc (idm,jdm))
      allocate(meanssh (idm,jdm))
      allocate(h     (idm,jdm))
      allocate(ficem (idm,jdm))
      allocate(hicem (idm,jdm))
      allocate(fice  (idm,jdm))
      allocate(hice  (idm,jdm))
      allocate(util  (idm*jdm+14))
      allocate(tmp   (idm,jdm))
      allocate(ut    (idm,jdm))
      allocate(vt    (idm,jdm))
      allocate(dp    (idm,jdm))
      allocate(tem   (idm,jdm))
      allocate(sal   (idm,jdm))
      allocate(speed (idm,jdm))
      allocate(mke   (idm,jdm))
      allocate(eke   (idm,jdm))
      allocate(ub    (idm,jdm))
      allocate(vb    (idm,jdm))
      allocate(dpsum (idm,jdm))

      if (idm>999 .or. jdm>999) then
         write(tag7,'(i5.5,a,i5.5)')idm,'x',jdm
      else
         write(tag7,'(i3.3,a,i3.3)')idm,'x',jdm
      end if

    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get file name, open them and  assign file units
      ficem=0.
      hicem=0.
      goticestate=.false.
      outname=''
      ub=0.
      vb=0.
      dpsum=0.
      open(10,file='multiproc.in')
      do i=1,maxfiles


         ! Check existence...
         read(10,*,end=100,err=100)filenames(i)
         inquire(file=trim(filenames(i)),exist=ex)

         ! Check it is a ".a"-file (wont work if name is changed....
         findab=index(filenames(i),'.a' ,back=.true.)

         if (.not.ex) then
            print *,'the file ',trim(filenames(i)),' does not exist'
            stop
         elseif (findab<=0) then
            ! You should only supply ".a" files in multiproc.in
            print *,'For new averaging file type -- only specify ".a" 
     &               type files in multiproc.in'
            stop
         endif


         ! Open file
         fileunits(i)=100+i
         call zaiopf(trim(filenames(i)),'old',fileunits(i))

         ! Open ".b" file as well
         fileb=trim(filenames(i))
         fileb(findab:findab+1)='.b'
         open(unit=fileunits(i),file=trim(fileb),status='old',
     &        form='formatted')

         print *,i,trim(filenames(i)),fileunits(i)

         if (i==1) then
            rungen=fileb(1:3)
         else if (rungen/=fileb(1:3)) then
            print *,'Warning: Files have different "rungen"'
         end if

      enddo
100   nrfiles=i-1
      close(10)
      if (nrfiles <= 0) stop 'multiproc.in is empty'
      if (nrfiles ==1) then
         find=index(filenames(1),'.',back=.true.)
         outname(1:find-1)=filenames(1)(1:find-1)
      end if
      if (nrfiles > 1) then
         outname(1:7)=filenames(1)(1:7)
         outname(8:10)='ave'
      endif


      ! Read header info...
      do i=1,nrfiles
       
         do j=1,12
            read(fileunits(i),*)
         end do

         ! Next entry is ave counter
         read(fileunits(i),*) avecounter(i)
         read(fileunits(i),*)

      end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main loop  - cycle files - read fields of each file. Make sure
! AVE-files contain excactly same number of arguments!
      utlevel=-1; vtlevel=-1; sallevel=-1; temlevel=-1;
      splevel=-1; dplevel=-1

      ficelevel=-1; hicelevel=-1
      tmp=0.
      ios=0
      maxlevel=-1
      inquire(iolength=j)fieldhead,util
      !write(*,*)'record length=',j
      open(10,file=trim(outname),status='replace',form='unformatted',
     &        access='direct',recl=j)
      print *,'Processing files now.. This may take a while if you have'
     &        //' many files'
      do while (ios==0)

         h_acc=0.
         counter_acc=0
         
         do i=1,nrfiles


            ! Read field info from .b file
            read (fileunits(i),117,iostat=ios) 
     &         c8,nstep,time,klevel,coord,hmin2,hmax2

            !if (i==1) print *,c8,klevel

            ! 
            if (ios/=0) goto 999

            ! For first pass, make sure we keep field name
            if (i==1) then 
               actual_field=c8
               actual_level=klevel
            else
               if (trim(actual_field)/=trim(c8) .or.
     &             klevel/=actual_level) then
                  print *,'Error: Fields or levels do not match '//
     &                    'in ave file'
                  print *,'First ave file has field:',trim(actual_field)
                  print *,'This  ave file has field:',trim(c8)
                  print *,'First ave file field at level:',actual_level
                  print *,'This  ave file field at level:',klevel
               end if
            end if


            
            ! Read field - bails out on error ...
            call zaiord(h,ip,.false.,hmin,hmax,fileunits(i))

            ! Accumulated field
            h_acc=h_acc+h
            counter_acc=counter_acc+avecounter(i)
            
         end do


         ! This routine converts ave-file ids to ids used in
         ! pakk-files... It also checks if field is 2d or 3d, and
         ! puts it int 2d and 3d var arrays
         print '(2a,i4)','actual_field,k:', actual_field,klevel
         call pakk_field_ids(trim(actual_field),klevel,pakk_id,
     &      char2d,char3d,n2d,n3d,maxfields)


         ! Dump data if pakk_id is ok  (known by routine pakk_field_ids)
         if (trim(pakk_id)/='EMPTY') then


            ! If one of 3d ocean vars, keep them here
            if (trim(pakk_id)=='TEM') then
               tem=h_acc
               temlevel=klevel
            elseif (trim(pakk_id)=='SAL') then
               sal=h_acc
               sallevel=klevel
            elseif (trim(pakk_id)=='UT') then
               ut=h_acc
               utlevel=klevel
            elseif (trim(pakk_id)=='VT') then
               vt=h_acc
               vtlevel=klevel
            elseif (trim(pakk_id)=='DP') then
               dp=h_acc
               dplevel=klevel
            elseif (trim(pakk_id)=='SPEED') then
               speed=h_acc
               splevel=klevel



            ! The following vars are assumed not to be
            ! layer-weighted...
            else

               klevel=max(klevel,1)
               call pakmsk(ip,h_acc/counter_acc,tmp,idm,jdm,util,length)
               call mkfldh(pakk_id,1,1,klevel,idm,idm,jdm,length,
     &                     fieldhead)
               write(10,rec=irec)fieldhead,util ; irec=irec+1

               ! Special case - dump meanssh
               if (trim(pakk_id)=='SSH') then
                  meanssh=0.
C                 print *,counter_acc
C                 print *,'min/max h_acc',
C    &               minval(h_acc),maxval(h_acc)
C                 meanssh=h_acc/counter_acc
                  do j=1,jdm
                  do i=1,idm
                     meanssh(i,j)=h_acc(i,j)/counter_acc
C                    print *,meanssh(i,j)
                  end do
                  end do
C                 print *,'min/max meanssh',
C    &               minval(meanssh),maxval(meanssh)
                  open(11,file='meanssh'//trim(tag7)//'.uf', 
     &                 status='unknown',form='unformatted')
                  write(11)meanssh
                  close(11)
C                 print *,'min/max meanssh',
C    &               minval(meanssh),maxval(meanssh)
               end if


            end if

            ! ICESTATE vars - produce average ice thickness /conc
            if (trim(pakk_id)=='FICE') then
               fice=h_acc/counter_acc
               ficelevel=klevel
               goticestate=.true.
            else if (trim(pakk_id)=='HICE') then
               hice=h_acc/counter_acc
               hicelevel=klevel
               goticestate=.true.
            end if


            ! Dump 3D ocean thickness -weighted vars when all layer vars
            ! are read
            if (klevel==temlevel .and. sallevel==klevel .and.
     &         utlevel==klevel   .and.  vtlevel==klevel .and.
     &         dplevel==klevel   .and.  splevel==klevel ) then

               tem=tem/(dp+.01)
               sal=sal/(dp+.01)
               ut =ut /(dp+.01)
               vt =vt /(dp+.01)
               speed =speed /(dp+.01)
               dp =dp / counter_acc

               ! Depth integrated velocity
               ub=ub+ut*dp
               vb=vb+vt*dp
               dpsum=dpsum+dp


               ! Get highest number of level of ocean model
               maxlevel=max(klevel,maxlevel)

               call pakmsk(ip,tem,tmp,idm,jdm,util,length)
               call mkfldh('TEM',1,1,klevel,idm,idm,jdm,length,
     &                     fieldhead)
               write(10,rec=irec)fieldhead,util ; irec=irec+1

               call pakmsk(ip,sal,tmp,idm,jdm,util,length)
               call mkfldh('SAL',1,1,klevel,idm,idm,jdm,length,
     &                     fieldhead)
               write(10,rec=irec)fieldhead,util ; irec=irec+1

               call pakmsk(ip,ut,tmp,idm,jdm,util,length)
               call mkfldh('UT',1,1,klevel,idm,idm,jdm,length,
     &                     fieldhead)
               write(10,rec=irec)fieldhead,util ; irec=irec+1

               call pakmsk(ip,vt,tmp,idm,jdm,util,length)
               call mkfldh('VT',1,1,klevel,idm,idm,jdm,length,
     &                     fieldhead)
               write(10,rec=irec)fieldhead,util ; irec=irec+1

               call pakmsk(ip,speed,tmp,idm,jdm,util,length)
               call mkfldh('SPEED',1,1,klevel,idm,idm,jdm,length,
     &                     fieldhead)
               write(10,rec=irec)fieldhead,util ; irec=irec+1

               call pakmsk(ip,dp,tmp,idm,jdm,util,length)
               call mkfldh('DP',1,1,klevel,idm,idm,jdm,length,
     &                     fieldhead)
               write(10,rec=irec)fieldhead,util ; irec=irec+1

               ! Mean and Eddy Kinetic energy
               mke=1.*0.5*(ut**2+vt**2)
               eke=1.*0.5*(speed-mke)



               ! EKE and MKE are derived fields
               call pakk_field_ids('EKE',klevel,pakk_id,
     &            char2d,char3d,n2d,n3d,maxfields)
               call pakk_field_ids('MKE',klevel,pakk_id,
     &            char2d,char3d,n2d,n3d,maxfields)

               call pakmsk(ip,eke,tmp,idm,jdm,util,length)
               call mkfldh('EKE',1,1,klevel,idm,idm,jdm,length,
     &                     fieldhead)
               write(10,rec=irec)fieldhead,util ; irec=irec+1

               call pakmsk(ip,mke,tmp,idm,jdm,util,length)
               call mkfldh('MKE',1,1,klevel,idm,idm,jdm,length,
     &                     fieldhead)
               write(10,rec=irec)fieldhead,util ; irec=irec+1

               utlevel=-1; vtlevel=-1; sallevel=-1; temlevel=-1;
               splevel=-1; dplevel=-1

               if (klevel==15) then
                  open(11,file='avetst.txt',form='formatted')
                  do j=1,jdm
                  do i=1,idm
                     write(11,'(2i5,6e14.4)') i,j, 
     &                  plon(i,j), plat(i,j), tem(i,j), sal(i,j),
     &                  dp  (i,j), h   (i,j)
                  end do
                  end do
                  close(11)
               end if
            elseif  (ficelevel==klevel.and.hicelevel==klevel) then

               !ICESTATE averaging 
               ficem=ficem+fice
               hicem=hicem+fice*hice
               ficelevel=-1
               hicelevel=-1
            end if

         end if
 999     continue

      end do



      ub=ub/(dpsum+.01)
      vb=vb/(dpsum+.01)
      klevel=0
      call pakk_field_ids('UBAVG',klevel,pakk_id,
     &      char2d,char3d,n2d,n3d,maxfields)
      call pakmsk(ip,ub,tmp,idm,jdm,util,length)
      call mkfldh(pakk_id,1,1,1,idm,idm,jdm,length,
     &                fieldhead)
      write(10,rec=irec)fieldhead,util ; irec=irec+1

      call pakk_field_ids('VBAVG',klevel,pakk_id,
     &      char2d,char3d,n2d,n3d,maxfields)
      call pakmsk(ip,vb,tmp,idm,jdm,util,length)
      call mkfldh(pakk_id,1,1,1,idm,idm,jdm,length,
     &                fieldhead)
      write(10,rec=irec)fieldhead,util ; irec=irec+1


      if (goticestate) then

         klevel=0
         call pakk_field_ids('ficem',klevel,pakk_id,
     &         char2d,char3d,n2d,n3d,maxfields)
         call pakmsk(ip,ficem,tmp,idm,jdm,util,length)
         call mkfldh(pakk_id,1,1,1,idm,idm,jdm,length,
     &                   fieldhead)
         write(10,rec=irec)fieldhead,util ; irec=irec+1

         call pakk_field_ids('hicem',klevel,pakk_id,
     &         char2d,char3d,n2d,n3d,maxfields)
         call pakmsk(ip,hicem,tmp,idm,jdm,util,length)
         call mkfldh(pakk_id,1,1,1,idm,idm,jdm,length,
     &                   fieldhead)
         write(10,rec=irec)fieldhead,util ; irec=irec+1
      end if

      close(10)




 117  format (a8,' = ',i11,f11.2,i3,f7.3,1p2e16.7)


      ! Close ave files
      do i=1,nrfiles
         close(fileunits(i))
         call zaiocl(fileunits(i))
      end do




      ! Dump file header
      varlist(    1:    n2d)=char2d(1:n2d)
      varlist(n2d+1:n2d+n3d)=char3d(1:n3d)

      if(mnproc==1) then
         write(lp,*)'Number of 2D fields is',n2d
         write(lp,*)'Number of 3D fields is',n3d
      end if

      print *
      print *,'2D Fields:'
      write(lp,'(10(a5,1x))')char2d(1:n2d)
      print *
      print *,'3D Fields:'
      write(lp,'(10(a5,1x))')char3d(1:n3d)
      !print *
      !write(lp,'(10(a5,1x))')varlist(1:n2d+n3d)


      call prepak26(rungen,outname,n2d,n3d,varlist,
     &             maxfields,idat,rt,sigma(1:maxlevel),idm,jdm,maxlevel)

      print *
      print *,'Average fields dumped to '//trim(outname)
      end program multiproc
