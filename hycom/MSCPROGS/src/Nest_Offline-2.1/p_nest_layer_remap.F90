      program nest_layer_remap
      !use m_layer_remapV1
      !use m_layer_remapV2
      use mod_sigma
      use m_layer_remapV4
      use m_layer_mixV1
      use m_dp0kini
      use m_parse_blkdat
      use m_get_nest_record
      use m_get_nest_info
      use m_read_nest_header
      use m_save_nest_header
      implicit none

      character(len= *),parameter :: fnestloc='Nest/nestloc.uf'
      character(len=80)           :: fndepths,nestingfile,nestingfileout
      character(len= 7)           :: tag7
      character(len=20)           :: ctitle

      ! Dimensions of inner grid
      integer :: idm, jdm, inest, &
         iidm, jjdm

      ! Arrays holding data on inner grid (where we need nesting
      ! conditions)
      real*8, allocatable, dimension(:,:) :: &
         inner_ior8
      real*4, allocatable, dimension(:,:) :: inner_ior4
      real, allocatable, dimension(:,:) :: tmpfld
      real, allocatable, dimension(:,:) :: &
         inner_lon, inner_lat, inner_depths
      integer, allocatable, dimension(:,:) :: &
         inner_ipiv, inner_jpiv

      real, allocatable, dimension(:,:,:) :: &
         old_int, new_int, &
         old_tem, new_tem, &
         old_sal, new_sal, &
         old_u  , new_u  , &
         old_v  , new_v  , &
         old_dens
      real, allocatable, dimension(:) :: targetdens, dp0k, targetint

      real, parameter  :: onem=9806.
      real, parameter  :: onemm=9806./1000.



         

      integer :: irec,klevel, inum_offset,ivar, &
         kdm, intvar, nrec_offset_in, num_offset_in, firstrec, lastrec,  &
         ibnd, ios, kdm_remap,lenstr
      real :: realvar
      logical :: ex, tstpoint, lfatal

      integer :: i1dim,i2dim,j1dim,j2dim
      integer :: i1,i2,j1,j2,indbnd,imax,day,year,hour,i,j,k,l,k2
      character(len=2) :: cbnd
      character(len=3) :: cvar

      integer :: kdminp, offsetinp, klevelinp, yearinp, dayinp, hourinp
      real :: maxvinp, minvinp, minvout, maxvout
      integer :: irec_out, ioffset,nrecout
      integer :: findcat, klist
      integer :: thflag

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
      allocate(inner_ior4(idm,jdm))
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

      fndepths='Nest/ndepths'//trim(tag7)//'.uf'

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

      ! Get target densities from blkdat.input (for nested model)
      call parse_blkdat('kdm   ','integer',realvar,kdm_remap,'blkdat.input')
      print *,'kdm_remap:',kdm_remap
      allocate(targetdens(kdm_remap))
      allocate(targetint (kdm_remap))
      allocate(dp0k      (kdm_remap))
      do k=1,kdm_remap
         call parse_blkdat('sigma ','real',targetdens(k),intvar,'blkdat.input',k)
         !print *,'target density:',k,targetdens(k)
      end do


      ! Get interface setup 
      call dp0kini(dp0k,kdm_remap,'blkdat.input')
      !print *,dp0k
      !stop


      ! Get some things from blkdat
      call parse_blkdat('thflag','integer',realvar,thflag,'blkdat.input')

      

      inquire(exist=ex, file='nestremap.in')
      if (.not.ex) then
         print *,'nestremap.in not found'
         stop '(nest_layer_remap)'
      end if
      open (89,file='nestremap.in',status='old')
      read(89,'(a80)',iostat=ios) nestingfile
      do while(ios==0)


         ! Check that only "i1" files are listed
         lenstr=len_trim(nestingfile)
         if (nestingfile(lenstr-1:lenstr)/='i1') then
            print *,' Only list "i1" files in nestremap.in'
            print *,'nestfile :',nestingfile
            stop '(nest_layer_remap)'
         else
            nestingfile=nestingfile(1:lenstr-3)
            !print *,nestingfile
         end if

         ! Check that the files are in a catalog named "Nest/
         findcat=index(nestingfile,'Nest/')
         if (findcat<=0) then
            print *,'Nesting files must be located in catalogue "Nest/"'
            stop '(nest_layer_remap)'
         end if

         ! Get kdm and offsets from header file
         call get_nest_info(trim(nestingfile),kdm,nrec_offset_in,num_offset_in)
         print *,'From get_nest_info:'
         print *,'kdm(number of layers)          :',kdm
         print *,'nrec_offset(records per time  ):',nrec_offset_in
         print *,'num_offset (time dumps in file):',num_offset_in


         ! For each boundary
         do ibnd=1,4
            
         if (ibnd==1) then 
            cbnd='i1'
            i1dim=1 ; i2dim=inest
            j1dim=1 ; j2dim=jdm
         elseif (ibnd==2) then 
            cbnd='ii'
            i1dim=iidm ; i2dim=idm
            j1dim=1 ; j2dim=jdm
         elseif (ibnd==3) then 
            cbnd='j1'
            i1dim=1 ; i2dim=idm
            j1dim=1 ; j2dim=inest
         elseif (ibnd==4) then 
            cbnd='jj'
            i1dim=1 ; i2dim=idm
            j1dim=jjdm ; j2dim=jdm
         end if

         ! Open old file
         inquire(exist=ex,file=trim(nestingfile)//'_'//cbnd)
         if (.not.ex) then
            print *,'No file '//trim(nestingfile)//'_'//cbnd
            stop '(nest_layer_remap)'
         end if
         

         write(6,'(a,a,a)') 'Reading records ',  &
            ' from file=',trim(nestingfile)//'_'//cbnd
         inquire(iolength=j)inner_ior4(i1dim:i2dim,j1dim:j2dim)
         open(10,file=trim(nestingfile)//'_'//cbnd,form='unformatted', &
              access='direct',recl=j,status='old')

         nestingfileout='Nest_remap/'//nestingfile(findcat+5:len_trim(nestingfile))
         write(6,'(a,a,a)') 'Dumping records ',  &
            ' to file=',trim(nestingfileout)//'_'//cbnd
         inquire(iolength=j)inner_ior4(i1dim:i2dim,j1dim:j2dim)
         open(20,file=trim(nestingfileout)//'_'//cbnd,form='unformatted', &
              access='direct',recl=j,status='replace',iostat=ios)
         if (ios/=0) then
            print *,'An error occured .. Did you create the catalogue '
            print *,'Nest_remap  before running this routine ?'
            stop '(nest_layer_remap)'
         end if

         ! Read records for all depth levels
         allocate(old_int(i1dim:i2dim,j1dim:j2dim,kdm))
         allocate(old_tem(i1dim:i2dim,j1dim:j2dim,kdm))
         allocate(old_sal(i1dim:i2dim,j1dim:j2dim,kdm))
         allocate(old_dens(i1dim:i2dim,j1dim:j2dim,kdm))
         allocate(old_u  (i1dim:i2dim,j1dim:j2dim,kdm))
         allocate(old_v  (i1dim:i2dim,j1dim:j2dim,kdm))

         allocate(new_int(i1dim:i2dim,j1dim:j2dim,kdm_remap))
         allocate(new_tem(i1dim:i2dim,j1dim:j2dim,kdm_remap))
         allocate(new_sal(i1dim:i2dim,j1dim:j2dim,kdm_remap))
         allocate(new_u  (i1dim:i2dim,j1dim:j2dim,kdm_remap))
         allocate(new_v  (i1dim:i2dim,j1dim:j2dim,kdm_remap))


         ! For each offset
         do ioffset=1,num_offset_in

            firstrec=(ioffset-1)*nrec_offset_in+1
            lastrec =(ioffset  )*nrec_offset_in

            ! Get depth interfaces
            do k=1,kdm
               call get_nest_record('INT',k,nestingfile,irec,firstrec,lastrec)
               read(10,rec=irec) inner_ior4(i1dim:i2dim,j1dim:j2dim)
               old_int(:,:,k)=inner_ior4(i1dim:i2dim,j1dim:j2dim)
            end do
            !print *,'ok'

            ! Get temp
            do k=1,kdm
               call get_nest_record('TEM',k,nestingfile,irec,firstrec,lastrec)
               read(10,rec=irec) inner_ior4(i1dim:i2dim,j1dim:j2dim)
               old_tem(:,:,k)=inner_ior4(i1dim:i2dim,j1dim:j2dim)
            end do
            !print *,'ok'

            ! Get temp
            do k=1,kdm
               call get_nest_record('SAL',k,nestingfile,irec,firstrec,lastrec)
               read(10,rec=irec) inner_ior4(i1dim:i2dim,j1dim:j2dim)
               old_sal(:,:,k)=inner_ior4(i1dim:i2dim,j1dim:j2dim)
            end do
            !print *,'ok'

            ! Get temp
            do k=1,kdm
               call get_nest_record('UT ',k,nestingfile,irec,firstrec,lastrec)
               read(10,rec=irec) inner_ior4(i1dim:i2dim,j1dim:j2dim)
               old_u  (:,:,k)=inner_ior4(i1dim:i2dim,j1dim:j2dim)
               !print *,'UT INPUT:',k,minval(old_u(:,:,k)),maxval(old_u(:,:,k))
            end do
            !print *,'ok'

            ! Get temp
            do k=1,kdm
               call get_nest_record('VT ',k,nestingfile,irec,firstrec,lastrec)
               read(10,rec=irec) inner_ior4(i1dim:i2dim,j1dim:j2dim)
               old_v  (:,:,k)=inner_ior4(i1dim:i2dim,j1dim:j2dim)
               !print *,'VT INPUT:',k,minval(old_v(:,:,k)),maxval(old_v(:,:,k))
            end do
            !print *,'ok'

            ! Deduce new interface values 
            do j=j1dim,j2dim
            do i=i1dim,i2dim


               ! --- Densities on old grid
               do k=1,kdm
                  if (thflag==0) then
                     old_dens(i,j,k)=sig0(old_tem(i,j,k),old_sal(i,j,k))
                  elseif (thflag==2) then
                     old_dens(i,j,k)=sig2(old_tem(i,j,k),old_sal(i,j,k))
                  else
                     print *,'Unknown thflag ',thflag
                     call exit(1)
                  end if
               end do

               !do k=1,kdm
               !  if (k>1) then
               !     if (old_dens(i,j,k-1)>old_dens(i,j,k)) then
               !        print *,'dens error:'
               !        print *,'old layer number:',i,j,k
               !        print '(a,i4,4f10.3)','tem, sal, dens,int ',k-1, &
               !                 old_tem(i,j,k-1),old_sal(i,j,k-1), &
               !                 old_dens(i,j,k-1),old_int(i,j,k-1)/onem
               !        print '(a,i4,4f10.3)','tem, sal, dens,int ',k, &
               !                 old_tem(i,j,k),old_sal(i,j,k), &
               !                 old_dens(i,j,k),old_int(i,j,max(k,1))/onem
               !        !print *
               !        !print *,'old temp:', old_tem(i,j,:)
               !        !print *,'old saln:', old_sal(i,j,:)
               !        !print *,'old dens:', old_dens(i,j,:)
               !        !print *,'old intf:', old_int(i,j,:)/onem
               !        !print *
               !        !stop '(nest_layer_remap:old density step)'
               !     end if
               !  end if
               !end do
               !print *,'test:',old_dens
                     

               !KAL -- calculate lowest mass-filled layer
               do k=2,kdm
                  if (abs(old_int(i,j,k)-old_int(i,j,k-1)) < onem ) then
                     ! Very simple
                     old_dens(i,j,k)=old_dens(i,j,k-1)
                  end if
               end do
               !print *,old_int(i,j,:)
               !print *,old_tem(i,j,:)
               !print *,old_sal(i,j,:)

                 

               ! set up target interfaces
!KAL           call layer_remapV1(old_int(i,j,:),old_tem(i,j,:),old_sal(i,j,:),kdm, &
!KAL                            targetdens,targetint,dp0k,kdm_remap)
!KAL -- This is the version used in the mercator remapping, seems to work better
               tstpoint=.false.
               !if (i==5 .and. j==100) tstpoint=.true.
               call layer_remapV4(old_int(i,j,:),old_dens(i,j,:),kdm, &
                                targetdens,targetint,dp0k,kdm_remap,thflag, &
                                inner_depths(i,j)*onem,tstpoint,lfatal)
               
               new_int(i,j,:)=targetint

               !if (i==5 .and. j==100) then
               !   print *,'-------'
               !   print *,old_int(i,j,:)
               !   print *
               !   print *,old_tem(i,j,:)
               !   print *
               !   print *,old_sal(i,j,:)
               !   print *
               !   print *,targetint
               !   print *,'-------'
               !end if



               ! Mix over new layers - temperature
               call layer_mixV1(old_int(i,j,:),old_tem(i,j,:),kdm, &
                                new_int(i,j,:),new_tem(i,j,:),kdm_remap)

               ! Mix over new layers - salinity
               call layer_mixV1(old_int(i,j,:),old_sal(i,j,:),kdm, &
                                new_int(i,j,:),new_sal(i,j,:),kdm_remap)

               ! NB - mixing does not conserve kinetic energy - it
               ! adheres to baroclinic considerations though


               ! Mix over new layers - u
               call layer_mixV1(old_int(i,j,:),old_u  (i,j,:),kdm, &
                                new_int(i,j,:),new_u  (i,j,:),kdm_remap)

               ! Mix over new layers - v
               call layer_mixV1(old_int(i,j,:),old_v  (i,j,:),kdm, &
                                new_int(i,j,:),new_v  (i,j,:),kdm_remap)

            end do
            end do

            !do k=1,kdm_remap
            !   print *,'UT remixed:',k,minval(new_u(:,:,k)),maxval(new_u(:,:,k))
            !end do
            !do k=1,kdm_remap
            !   print *,'VT remixed:',k,minval(new_v(:,:,k)),maxval(new_v(:,:,k))
            !end do


            ! Almost done Now we just have to dump the data

            print *,' |--> Treating time offset ',ioffset
            if (ioffset==1) irec_out=0 ! Only on first "offset" pass
            do irec=firstrec,lastrec


               ! Read header of input record
               call read_nest_header(nestingfile,kdminp,offsetinp,cvar,klevelinp, &
                                     yearinp,dayinp,hourinp,minvinp,maxvinp,irec)

               if (offsetinp/=ioffset) then 
                  print *,'Offset error'
                  stop '(nest_layer_remap)'
               end if

               ! Dump interface values to OUTPUT file if we reached first record
               ! of interface values
               if (trim(cvar)=='INT' .and. klevelinp==1) then
                  !print *,'INT dumping starts at output record ',irec_out+1
                  do k=1,kdm_remap
                     irec_out=irec_out+1
                     inner_ior4(i1dim:i2dim,j1dim:j2dim)=new_int(:,:,k)
                     minvout=minval(inner_ior4(i1dim:i2dim,j1dim:j2dim))
                     maxvout=maxval(inner_ior4(i1dim:i2dim,j1dim:j2dim))
                     write(20,rec=irec_out) inner_ior4(i1dim:i2dim,j1dim:j2dim)
                     if (ibnd==1) then !only on second boundary pass...
                        call save_nest_header(nestingfileout,kdm_remap,offsetinp,cvar,k, &
                                              yearinp,dayinp,hourinp,minvout,maxvout,irec_out)
                     end if
                     !print *,irec,irec_out,'INT'
                  end do
               else if (trim(cvar)=='INT' .and. klevelinp>1) then
                  cycle


               ! Dump interface values to OUTPUT file if we reached first record
               ! of temperature values
               else if (trim(cvar)=='TEM' .and. klevelinp==1) then
                  !print *,'TEM dumping starts at output record ',irec_out+1
                  do k=1,kdm_remap
                     irec_out=irec_out+1
                     inner_ior4(i1dim:i2dim,j1dim:j2dim)=new_tem(:,:,k)
                     minvout=minval(inner_ior4(i1dim:i2dim,j1dim:j2dim))
                     maxvout=maxval(inner_ior4(i1dim:i2dim,j1dim:j2dim))
                     write(20,rec=irec_out) inner_ior4(i1dim:i2dim,j1dim:j2dim)
                     if (ibnd==2) then !only on second boundary pass...
                        call save_nest_header(nestingfileout,kdm_remap,offsetinp,cvar, &
                                              k,yearinp,dayinp,hourinp,minvout,maxvout,irec_out)
                     end if
                     !print *,irec,irec_out,'TEM'
                  end do
               else if (trim(cvar)=='TEM' .and. klevelinp>1) then
                  cycle
               

               ! Dump interface values to OUTPUT file if we reached first record
               ! of salinity values
               else if (trim(cvar)=='SAL' .and. klevelinp==1) then
                  !print *,'SAL dumping starts at output record ',irec_out+1
                  do k=1,kdm_remap
                     irec_out=irec_out+1
                     inner_ior4(i1dim:i2dim,j1dim:j2dim)=new_sal(:,:,k)
                     minvout=minval(inner_ior4(i1dim:i2dim,j1dim:j2dim))
                     maxvout=maxval(inner_ior4(i1dim:i2dim,j1dim:j2dim))
                     write(20,rec=irec_out) inner_ior4(i1dim:i2dim,j1dim:j2dim)
                     if (ibnd==2) then !only on second boundary pass...
                        call save_nest_header(nestingfileout,kdm_remap,offsetinp,cvar, &
                                              k,yearinp,dayinp,hourinp,minvout,maxvout,irec_out)
                     end if
                     !print *,irec,irec_out,'SAL',maxval(abs(inner_ior4(i1dim:i2dim,j1dim:j2dim)))
                  end do
               else if (trim(cvar)=='SAL' .and. klevelinp>1) then
                  cycle
              

               ! Dump interface values to OUTPUT file if we reached first record
               ! of baroclinic u values
               elseif (trim(cvar)=='UT' .and. klevelinp==1) then
                  !print *,'UT  dumping starts at output record ',irec_out+1
                  do k=1,kdm_remap
                     irec_out=irec_out+1
                     inner_ior4(i1dim:i2dim,j1dim:j2dim)=new_u  (:,:,k)
                     minvout=minval(inner_ior4(i1dim:i2dim,j1dim:j2dim))
                     maxvout=maxval(inner_ior4(i1dim:i2dim,j1dim:j2dim))
                     write(20,rec=irec_out) inner_ior4(i1dim:i2dim,j1dim:j2dim)
                     if (ibnd==2) then !only on second boundary pass...
                        call save_nest_header(nestingfileout,kdm_remap,offsetinp,cvar,k,yearinp,&
                                              dayinp,hourinp,minvout,maxvout,irec_out)
                     end if
                     !print *,irec,irec_out,'UT',maxval(abs(inner_ior4(i1dim:i2dim,j1dim:j2dim)))
                  end do
               else if (trim(cvar)=='UT' .and. klevelinp>1) then
                  cycle

               ! Dump interface values to OUTPUT file if we reached first record
               ! of baroclinic v values
               else if (trim(cvar)=='VT' .and. klevelinp==1) then
                  !print *,'VT  dumping starts at output record ',irec_out+1
                  do k=1,kdm_remap
                     irec_out=irec_out+1
                     inner_ior4(i1dim:i2dim,j1dim:j2dim)=new_v  (:,:,k)
                     minvout=minval(inner_ior4(i1dim:i2dim,j1dim:j2dim))
                     maxvout=maxval(inner_ior4(i1dim:i2dim,j1dim:j2dim))
                     write(20,rec=irec_out) inner_ior4(i1dim:i2dim,j1dim:j2dim)
                     if (ibnd==2) then !only on second boundary pass...
                        call save_nest_header(nestingfileout,kdm_remap,ioffset,cvar,k,yearinp, &
                                              dayinp,hourinp,minvout,maxvout,irec_out)
                     end if
                     !print *,irec,irec_out,'VT',maxval(abs(inner_ior4(i1dim:i2dim,j1dim:j2dim)))
                  end do
               else if (trim(cvar)=='VT' .and. klevelinp>1) then
                  cycle


               ! All other conditionals simply duplicate input to output -- May
               ! need to add ecosys vars to the logic above
               else
                  !print *,cvar//' dumping starts at output record ',irec_out+1
                  irec_out=irec_out+1
                  read (10,rec=irec    ) inner_ior4(i1dim:i2dim,j1dim:j2dim)
                  !print *,irec,irec_out,cvar
                  minvout=minval(inner_ior4(i1dim:i2dim,j1dim:j2dim))
                  maxvout=maxval(inner_ior4(i1dim:i2dim,j1dim:j2dim))
                  write(20,rec=irec_out) inner_ior4(i1dim:i2dim,j1dim:j2dim)
                  if (ibnd==2) then !only on second boundary pass...
                     call save_nest_header(nestingfileout,kdm_remap,offsetinp,cvar,klevelinp,yearinp, &
                                           dayinp,hourinp,minvout,maxvout,irec_out)
                  end if
               end if

               !print *,irec,irec_out

            end do

            ! irec_out now is equal to num records in one time
            if (ioffset==1) nrecout=irec_out

         end do ! ioffset

         close(10) 
         close(20)

         deallocate(old_dens)
         deallocate(old_int)
         deallocate(old_tem)
         deallocate(old_sal)
         deallocate(old_u  )
         deallocate(old_v  )

         deallocate(new_int)
         deallocate(new_tem)
         deallocate(new_sal)
         deallocate(new_u  )
         deallocate(new_v  )

         end do ! ibnd
         print '(a,i4.4,a,i3.3)','Finito for year ',yearinp, ' and day ',dayinp
         print '(a,i2,a,i2,a)','Old file had ',kdm,' layers, new file has ',kdm_remap,' layers'
         print '(a,i4,a,i4,a)','Input  file had ',num_offset_in,' times and ',nrec_offset_in,' records per time'
         print '(a,i4,a,i4,a)','Output file has ',num_offset_in,' times and ',nrecout       ,' records per time'
         print *
         print *
         print *
         read(89,'(a80)',iostat=ios) nestingfile
      end do ! files
      close(89)




      end program nest_layer_remap
