      program gen_ens
#define ICE 
      ! Reads initial model field from hycom input files, and (if present)
      ! from ice input file (old format). New hycom files are created
      ! by perturbing the layer thickness and ice thickness
      use mod_xc
      use mod_za
      use mod_grid
      use mod_pseudo
      use m_parse_blkdat
      implicit none 

      character(len=80) infile,outfile,filebase,cline,rstfile,icefile
      character(len=6) cvarin
      character(len=7) tag7
      character(len=8) var
      character(len=3) cmem
      logical ex, lok, l_keep

      real rv,rh,alp,bet, sd_d
      integer iens,i,j,ensdim, n1,n2, layer, level
      integer :: kdm
      integer :: ios,fab

      real :: hmina,hmaxa,hminb,hmaxb,rvar, rdummy

      character(len=3) rungen, censdim
      character(len=50) :: tmparg
      integer k,l,iostat
      real, allocatable, dimension(:,:) :: fld2d,ranfld,accranfld, wrongtopo
      real*8, allocatable, dimension(:,:) :: iofld
      real, allocatable, dimension(:,:,:) :: dptmp
      integer, allocatable :: dummy(:,:)
      real*8, dimension(:,:), allocatable ::  ficem,hicem,hsnwm,ticem,tsrfm

      real, parameter :: onem=9806.


      inquire( exist = ex , file ='gen_ens.in') 

      if (ex) then
         open(11,file='gen_ens.in',action='read', status='old')
         read(11,*) rstfile
         read(11,*) ensdim ! # ensemble members
         read(11,*) l_keep !
         read(11,*) sd_d   ! Std dev of log(d), no unit, 0.1 = 10%
         read(11,*) rv     ! Vertical correlation range (in nb layers)
         read(11,*) rh     ! horizontal correlation range (in nb of grid cells)
         close(11)

      else

         ! print informative message
         print *,'This routine will create perturbations of the  hycom layer'
         print *,'thickness using spatial fields with prescribed correlation'
         print *,'in the horizontal and vertical. Correlation lengths are in'
         print *,'layers (vertical) and grid cells (horizontal).'
         print *

         print *,'Start file  (member 1), correlation scales and standard '
         print *,'deviation of perturbations are specified in the input file'
         print *,'gen_ens.in An example of this file is found  '
         print *,'between the ----- markers below. First line is the restart'
         print *,'file to use.'
         print *,'-----'
         print *,'FORrestart2008_280_00.a'
         print *,'50           # Total number of ensemble members to create'
         print *,'F            # just set to false..'
         print *,'0.05         # Std dev of log(d), no unit, 0.1 = 10%'
         print *,'3            # Vertical correlation range   (layers)'
         print *,'10           # Horizontal correlation range (grid cells)'
         print *,'-----'
         print *

         print *,'If a ICE restartfile is found, this will be perturbed as well, but '
         print *,'this time by perturbing the ice thickness fields. (In the '
         print *,'same way, if a field "hicem" is found in the hycom restart file,'
         print *,'then this field will be perturbed)'
         print *
         print *,'gen_ens.in missing'
         call exit(1)
      end if


      fab=max(index(rstfile,'.a'),index(rstfile,'.b'))
      filebase=rstfile(1:fab-1)

      ! Initialize
      call parse_blkdat('kdm   ','integer',rdummy,kdm)
      call xcspmd()
      call zaiost()
      call get_grid()
      allocate(fld2d    (idm,jdm))
      allocate(dummy    (idm,jdm))
      allocate(ranfld   (idm,jdm))
      allocate(accranfld(idm,jdm))
      allocate(wrongtopo(idm,jdm))
      allocate(iofld    (idm,jdm))
      allocate(dptmp    (idm,jdm,kdm))


      lp=6

      print '(a,3f10.2)','scales: ', sd_d, rv, rh


!      sd_d = 0.1  ! Standard deviation of log(d) no unit
!                  ! remember enstmp std Gaussian : Prob(|enstmp|>2)= 5% 
!      rv=2.0e00   ! sort of vertical correlation range (in nb of layers)
!      rh=25.00   ! horizontal decorrelation scale (number of grid cells)

      !  Generates the vertical correlation of the ensemble.
      alp=exp(-1.0/rv)
      print *,'Alp=',alp
      bet=sqrt(1.0-alp**2) ! keeps var(enstmp)=var(ensmem)
      print *,'Bet=',bet

      ! original Hycom state kept as first member in record 1 - 

      ! FFT dimensions
      n1=2**(ceiling(log(float(idm))/log(2.)))
      n2=2**(ceiling(log(float(jdm))/log(2.)))

      if(l_keep) then 
         print '(a)','Re-Using ensemble members - HYCOM:'
      else
         print '(a)','Creating new ensemble members - HYCOM:'
      endif

      do iens=2,ensdim  

         ! set file names
         write(cmem,'(i3.3)') iens
         if(l_keep) then 
            infile=trim(filebase)//'_mem'//cmem
            outfile=trim(filebase)//'_tmp'//cmem
         else
            infile=trim(filebase)
            outfile=trim(filebase)//'_mem'//cmem
         endif

         print *
         print '(a,i4,a,a)','Member ',iens, ' file          : ',trim(outfile)
         print '(a)','Fields will be based on   : '//infile

         wrongtopo=0.

         ! Open infile
         call zaiopf(trim(infile)//'.a', 'old', 701)
         open (unit=701,file=trim(infile)//'.b',status='old')

         ! Open outfile
         call zaiopf(trim(outfile)//'.a', 'replace', 801)
         open (unit=801,file=trim(outfile)//'.b',status='replace')

         ! Start reading infile - 2 first lines are header
         read(701,'(a80)') cline
         write(801,'(a80)') cline
         read(701,'(a80)') cline
         write(801,'(a80)') cline


         ! Start reading - until we reach "dp" fields, this is 
         ! just copying fields from input to output
         do while (.not.(trim(var)=='v'.and.level==2 .and. layer==kdm))
            read(701,'(a80)') cline
            var=cline(1:8)
            i = index(cline,'=')
            read (cline(i+1:),*) layer,level,hminb,hmaxb
            !print *,var,layer,level,hminb,hmaxb
            write(801,4100) var,layer,level,hminb,hmaxb
            call zaiord(fld2d, dummy,.false., hmina,hmaxa,  701)
            call zaiowr(fld2d, dummy,.false., hmina,hmaxa,  801, .false.)
         end do


         ! Process dp layers -- dp must be the same for both time steps,
         ! otherwise temporal gradients may cause problems
         print '(a)','-perturbing layers '
         do k=1,kdm 

            read(701,'(a80)') cline
            var=cline(1:8)
            i = index(cline,'=')
            read (cline(i+1:),*) layer,level,hminb,hmaxb
            !print *,var,layer,level,hminb,hmaxb
            call zaiord(fld2d, dummy,.false., hmina,hmaxa,  701)

            if (trim(var)=='dp' .and. level == 1 .and. layer==k) then


               ! 3D vertically correlated
               where(fld2d>1e28) fld2d=0.
               call pseudo2D(   ranfld,idm,jdm,1,rh,n1,n2)
               if (k==1) then
                  call pseudo2D(accranfld,idm,jdm,1,rh,n1,n2)
               else
                  accranfld=alp*accranfld+bet*ranfld
               end if
               !print *,'dp test unperturbed',k,minval(fld2d)/onem,maxval(fld2d)/onem
               fld2d=fld2d*exp(accranfld*sd_d-0.5*sd_d**2)
               !print *,'dp test perturbed  ',k,minval(fld2d)/onem,maxval(fld2d)/onem

               ! Accumulated depth -- needed later
               wrongtopo=wrongtopo+fld2d

               ! We can afford one 3D variable (about 60 MB for a 800x880x22 grid)
               dptmp(:,:,k)=fld2d
            end if
         end do
         print '(a,2f14.2)','topo test dpsum  perturbed:',minval(wrongtopo)/onem,maxval(wrongtopo)/onem


         ! layer sum may be wrong, this corrects it
         do k=1,kdm 
            do j=1,jdm
            do i=1,idm
            if (depths(i,j)>.1) then
               dptmp(i,j,k)=dptmp(i,j,k)*depths(i,j)*onem/wrongtopo(i,j)
            else
               dptmp(i,j,k)=0.
            end if
            end do
            end do
         end do
         print '(a,2f14.2)','topo test dpsum  corrected:',minval(sum(dptmp,3))/onem,maxval(sum(dptmp,3))/onem
         print '(a,2f14.2)','topo test depths          :',minval(depths),maxval(depths)


         ! Spit new layer thickness into outfile
         do k=1,kdm 
            call zaiowr(dptmp(:,:,k), dummy,.false., hmina,hmaxa,  801, .false.)
            write(801,4100) var,k,level,hmina,hmaxa
         end do


         ! Spit same layer thickness into outfile's second time level
         do k=1,kdm 
            ! Need to skip forward in infile as well
            read(701,'(a80)') cline
            var=cline(1:8)
            i = index(cline,'=')
            read (cline(i+1:),*) layer,level,hminb,hmaxb
            call zaiord(fld2d, dummy,.false., hmina,hmaxa,  701)
            call zaiowr(dptmp(:,:,k), dummy,.false., hmina,hmaxa,  801, .false.)
            write(801,4100) var,layer,level,hmina,hmaxa
         end do



         ! The rest of the outfile is the same as the infile  except for ice
         ! thickness (if present)
         ios = 0
         do while (ios==0)
            read(701,'(a80)',iostat=ios) cline
            if (ios==0) then
               var=cline(1:8)
               i = index(cline,'=')
               read (cline(i+1:),*) layer,level,hminb,hmaxb
               call zaiord(fld2d, dummy,.false., hmina,hmaxa,  701)

               ! If this is a ice thickness field, perturb that as well
               if (trim(var)=='hicem') then
                  print '(a)','-perturbing ice thickness '
                  call pseudo2D(   ranfld,idm,jdm,1,rh,n1,n2)
                  where(fld2d>1e28) fld2d=0.
                  print '(a,2f14.2)','old max ice thickness :',maxval(fld2d)
                  fld2d=fld2d*exp(ranfld*sd_d-0.5*sd_d**2) ! Same approach as for layer thickness
                  print '(a,2f14.2)','new max ice thickness :',maxval(fld2d)
               end if
               call zaiowr(fld2d, dummy,.false., hmina,hmaxa,  801, .false.)
               write(801,4100) var,layer,level,hmina,hmaxa
            end if
         end do

         ! Close up shop
         call zaiocl(701)
         close(701)
         call zaiocl(801)
         close(801)
         if(l_keep) then 
            print*, 'Overwrite tmp member'
            call system('mv '//trim(outfile)//'.a '//trim(infile)//'.a')
            call system('mv '//trim(outfile)//'.b '//trim(infile)//'.b')
         endif
      enddo

      
      deallocate(fld2d)
      deallocate(dummy)
      deallocate(accranfld)
      deallocate(depths)
      deallocate(wrongtopo)
      deallocate(iofld)
      deallocate(dptmp)

#ifdef ICE
! ICE ensemble
      print *
      print *
      icefile=trim(filebase)//'ICE.uf'
      inquire(exist=ex, file=trim(icefile))
      if (.not. ex) then 
         print *,'No ice file present'
         stop
      end if

      allocate(ficem(idm,jdm))
      allocate(hicem(idm,jdm))
      allocate(hsnwm(idm,jdm))
      allocate(ticem(idm,jdm))
      allocate(tsrfm(idm,jdm))

      inquire(iolength=j) ficem,hicem,hsnwm,ticem,tsrfm  
      open(10,file=trim(icefile),status='old',form='unformatted', &
              access='direct',recl=j)

      rh=25. ; sd_d=0.1 ! By authority, should not be the same as ocean fields
      print *,'Creating ensemble members - ICE: rh=25. ; sd_d=0.1'
      do iens=2,ensdim
!        No -keep for ICE, all members 2..100 are replaced.
         read(10 ,rec=1   ,iostat=ios) ficem,hicem,hsnwm,ticem,tsrfm
         call pseudo2D(ranfld,idm,jdm,1,rh,n1,n2)
         hicem=hicem*exp(ranfld*sd_d-0.5*sd_d**2)
         write(10,rec=iens,iostat=ios) ficem,hicem,hsnwm,ticem,tsrfm
      enddo

      deallocate(ficem)
      deallocate(hicem)
      deallocate(hsnwm)
      deallocate(ticem)
      deallocate(tsrfm)
#endif
      deallocate(ranfld)

      close(10)  ! blkdat.input

 4100 format(a,': layer,tlevel,range = ',i3,i3,2x,1p2e16.7)

      end program gen_ens


