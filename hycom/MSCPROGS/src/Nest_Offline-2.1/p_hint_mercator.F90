! KAL -- This routine interpolates horizontal mercator fields to a NERSC model
! KAL -- No vertical interpolation done by this routine.
! KAL -- 
! KAL -- Usable/Inital version: 31.01.2007 - Knut Lisæter



      program mercator_horizontal_interpolation
      use netcdf
      use mod_xc
      use mod_za
      use mod_grid, only: hyclat => plat, hyclon => plon, hycbathy => depths, &
                    get_grid, hyclonU => ulon, hyclatU => ulat, &
                    hyclonV => vlon, hyclatV => vlat
      use mod_mercatorgrid
      use mod_confmap
      use m_nearestpoint
      implicit none


       character(len=80) :: testfile, ctemp, csaln, cuvel, cvvel,  testfileu, testfileV,cssh, &
          filebase
       integer :: ierr
       integer :: k,ivar
       integer :: ncid, ncidU, ncidV
       integer :: tempid,zvarid,zvaridU,lonid,latid, tmpid, salnid, uvelid, vvelid, zvaridV, &
                  sshid, depthid

       real, allocatable :: temp2d(:,:),  &
           tmpr(:,:),  &
           tmprmask(:,:), temp2dold(:,:), tmprold(:,:), &
           hycnearest(:,:), ubavg(:,:),vbavg(:,:), ssh(:,:), tmpssh(:,:), &
           newdepth(:,:)
        real*8, allocatable :: iofld(:,:)
       integer, allocatable :: ip(:,:),maxk(:,:)
       integer, dimension(nf90_max_var_dims) :: dimIDs
       integer :: ipivtmp,jpivtmp,ipivh,jpivh
       integer :: i2,j2,i,j
       integer :: ndims
       real ::  amax, amin, dist, mindist
       real, allocatable :: zlev(:), zlevU(:), zlevV(:)
       logical :: usefile,ass,ex
       real :: a1,a2,a3,a4
       real :: dptmp
       real :: missing_value
       character(len=8)  :: chfld
       character(len=11) :: tag11
       real, dimension(:,:), allocatable :: mercmasktst

       real    :: tmpdnear, missing_value_ssh
       integer :: knear

       real :: reflength
       logical :: maskok, exT, exU, exV

       real, external :: spherdist
#if defined(IARGC)
       integer*4, external :: iargc
#endif


       if ( iargc() /= 1 ) then
          print *,'Usage:'
          print *,'   hint_mercator basename'
          print *
          print *,'base name is name of mercator file up to and including "grid"'
          print *,'Example -- basename equal to ext-NATL4-T17_y2005m01d01_grid '
          print *,'will use the following files:'
          print *,'             ext-NATL4-T17_y2005m01d01_gridT.nc'
          print *,'             ext-NATL4-T17_y2005m01d01_gridU.nc'
          print *,'             ext-NATL4-T17_y2005m01d01_gridV.nc'
          call exit(1)
       else
          call getarg(1,filebase)
       end if

       testfile =trim(filebase)//'T.nc';
       testfileU=trim(filebase)//'U.nc';
       testfileV=trim(filebase)//'V.nc';

       inquire(exist=exT,file=trim(testfile ))
       inquire(exist=exU,file=trim(testfileU))
       inquire(exist=exV,file=trim(testfileV))

       if (.not. (exT.and.exU.and.exV)) then
          print *,'One of these files are missing :'
          print *,trim(testfile )
          print *,trim(testfileU)
          print *,trim(testfileV)
          call exit(1)
       end if


       ! Read Mercator netcdf file 
       !testfile='ext-NATL4-T17_y2005m01d01_gridT.nc'
       !testfileU='ext-NATL4-T17_y2005m01d01_gridU.nc'
       !testfileV='ext-NATL4-T17_y2005m01d01_gridV.nc'
       ctemp='votemper'
       csaln='vosaline'
       cuvel='vozocrtx'
       cvvel='vomecrty'
       cssh ='sossheig'
       ierr=NF90_OPEN(testfile,NF90_NOWRITE,ncid)
       ierr=NF90_OPEN(testfileU,NF90_NOWRITE,ncidU)
       ierr=NF90_OPEN(testfileV,NF90_NOWRITE,ncidV)


       ! Get HYCOM lon lat and depths
       ! Initialize Arrai IO
       call xcspmd()
       call zaiost()
       call get_grid()
       allocate(ip(idm,jdm))
       print *,'hyc lon:',minval(hyclon),maxval(hyclon)
       print *,'hyc lat:',minval(hyclat),maxval(hyclat)


       ! This initializes the mercator grid size lengths
       call mercgrid_init(idm,jdm)
       call mercgrid_mask(1)

       inquire(exist=ex,file='regional.grid.a')
       if (.not.ex) stop '(p_hint_mercator: regional.grid.a not present)'
       inquire(exist=ex,file='regional.depth.a')
       if (.not.ex) stop '(p_hint_mercator: regional.depth.a not present)'


       ! Init confmap
       call initconfmap(idm,jdm)

       usefile=.true.

       call mercgrid_pivots(merclon,merclat,hyclon,hyclat, &
                            inear,jnear,distnear,idm,jdm,usefile,'t')

       call mercgrid_pivots(merclonu,merclatu,hyclonu,hyclatu, &
                            inearu,jnearu,distnearu,idm,jdm,usefile,'u')

       call mercgrid_pivots(merclonv,merclatv,hyclonv,hyclatv, &
                            inearv,jnearv,distnearv,idm,jdm,usefile,'v')


       reflength=maxval(distnear)/4
       print *,'reflength is ', reflength





       ! Get varid
       ierr=nf90_inq_varid(ncid,'deptht',zvarid)
       ierr=nf90_inq_varid(ncid,trim(ctemp),tempid)
       ierr=nf90_inq_varid(ncid,trim(csaln),salnid)
       ierr=nf90_inq_varid(ncid,trim(cssh) ,sshid)
       ierr=nf90_inq_varid(ncid,trim('deptht') ,depthid)
       write(6,'(a)') NF90_STRERROR(ierr)
       !ierr=nf90_inquire_variable(ncid,tempid,dimids=dimids,ndims=ndims)
       if (ierr/=NF90_NOERR) stop '(error getting v temp)'


       ! Get varid
       ierr=nf90_inq_varid(ncidU,'depthu',zvaridU)
       ierr=nf90_inq_varid(ncidU,trim(cuvel),uvelid)
       !ierr=nf90_inquire_variable(ncidU,uvelid,dimids=dimids,ndims=ndims)
       if (ierr/=NF90_NOERR) stop '(error getting v uvel)'
       
       ! Get varid
       ierr=nf90_inq_varid(ncidV,'depthv',zvaridV)
       ierr=nf90_inq_varid(ncidV,trim(cvvel),vvelid)
       !ierr=nf90_inquire_variable(ncidV,vvelid,dimids=dimids,ndims=ndims)
       if (ierr/=NF90_NOERR) stop '(error getting v uvel)'

       ! Allocate 2D vars to hold data
       allocate(tmpr   (idm,jdm))
       allocate(tmprold(idm,jdm))
       call zaiopf('merchint.a','replace',99)
       open (98,file='merchint.b',status='replace')

       tmpr=real(inear(:,:,1))
       call zaiowr(tmpr,ip,.false.,amin,amax,99,.false.) !plon
       write(98,104) 'ipiv    ',0,0.,amin,amax
       tmpr=real(jnear(:,:,1))
       call zaiowr(tmpr,ip,.false.,amin,amax,99,.false.) !plon
       write(98,104) 'jpiv    ',0,0.,amin,amax


       allocate(tmpssh(dimlenx,dimleny))
       allocate(ssh(idm,jdm))
       allocate(newdepth(idm,jdm))
       allocate(ubavg(idm,jdm))
       allocate(vbavg(idm,jdm))
       allocate(temp2d(dimlenx,dimleny))
       allocate(temp2dold(dimlenx,dimleny))
       allocate(mercmasktst(dimlenx,dimleny))
       allocate(tmprmask(idm,jdm))
       allocate(hycnearest(idm,jdm))
       allocate(maxk(idm,jdm))
       allocate(zlev (dimlenz))
       allocate(zlevU(dimlenz))
       allocate(zlevV(dimlenz))
       maxk=0
       ssh=0.
       newdepth=0.
       ubavg=0.
       vbavg=0.

       ! z level of model
       ierr=nf90_get_var(ncid,zvarid,zlev)

       ! Go through all Mercator depths 
       do ivar=1,4
       do k=1,dimlenz


          ! Use temperatures for testing
          ! Get "undef" values
          if (ivar==1) then
             ierr=nf90_get_var(ncid,tempid,temp2d,start=(/1,1,k,1/),count=(/dimlenx,dimleny,1,1/))
             ierr=nf90_get_att(ncid,tempid,'missing_value',missing_value)
             chfld='temp    '


             ! Get ssh and deptht from the same file 
             if (k==1) then
                ierr=nf90_get_var(ncid,sshid,tmpssh,start=(/1,1,1/),count=(/dimlenx,dimleny,1/))
                ierr=nf90_get_att(ncid,sshid,'missing_value',missing_value_ssh)
             end if

          elseif (ivar==2) then
             ierr=nf90_get_var(ncid,salnid,temp2d,start=(/1,1,k,1/),count=(/dimlenx,dimleny,1,1/))
             ierr=nf90_get_att(ncid,salnid,'missing_value',missing_value)
             chfld='saln    '
          elseif (ivar==3) then
             ierr=nf90_get_var(ncidu,uvelid,temp2d,start=(/1,1,k,1/),count=(/dimlenx,dimleny,1,1/))
             ierr=nf90_get_att(ncidu,uvelid,'missing_value',missing_value)
             chfld='utot    '
          elseif (ivar==4) then
             ierr=nf90_get_var(ncidV,vvelid,temp2d,start=(/1,1,k,1/),count=(/dimlenx,dimleny,1,1/))
             ierr=nf90_get_att(ncidV,vvelid,'missing_value',missing_value)
             chfld='vtot    '
          else
             print *,'ivar range ...'
             stop
          end if

          if (k==1) print *,'Processing '//chfld

          ! mask for t-cells (temperature)
          !call mercgrid('tmask',mercmask,k)
          !mercmask=1
          !where(temp2d==missing_value)mercmask=0

          ! Mask for each z-level
          call mercgrid_mask(k)
          mercmasktst=1
          where (temp2d==missing_value) mercmasktst=0

          !if (ivar==1.or.ivar==2) then
          !   print *,count(mercmasktst/=mercmask)
          !elseif (ivar==3) then
          !   print *,count(mercmasktst/=mercmasku)
          !elseif (ivar==4) then
          !   print *,count(mercmasktst/=mercmaskv)
          !end if

          ! For some reason we have to  do this...
          mercmask=mercmasktst



          !print *,zlev(k),minval(temp2d),maxval(temp2d)

          do j=1,jdm
          do i=1,idm
             i2=inear(i,j,1)
             j2=jnear(i,j,1)
             tmpr(i,j)=temp2d(i2,j2)
             tmprmask(i,j)=mercmask(i2,j2)

             !! sum (fld*e^(-x/max(x))) / sum( e^(-x/max(x)))
             tmpr(i,j) =  0.
             maskok=.true.
             do knear=1,numnear
                i2=inear(i,j,knear)
                j2=jnear(i,j,knear)
                tmpdnear=distnear(i,j,knear)
                tmpr(i,j)= tmpr(i,j) + temp2d(i2,j2) * exp( - tmpdnear / reflength )

                if (ivar==1.and.k==1) then
                   ssh(i,j)     =      ssh(i,j) + tmpssh(i2,j2) * exp( - tmpdnear / reflength )
                   newdepth(i,j)= newdepth(i,j) + mercbathy(i2,j2) * exp( - tmpdnear / reflength )
                end if



                maskok=maskok.and.mercmask(i2,j2)==1
             end do

             if (maskok) then
                tmpr(i,j)=tmpr(i,j) / sum( exp( - distnear(i,j,:) / reflength ) ) 
                if (ivar==1.and.k==1) then
                   ssh (i,j)= ssh(i,j) / sum( exp( - distnear(i,j,:) / reflength ) ) 
                   newdepth (i,j)= newdepth(i,j) / sum( exp( - distnear(i,j,:) / reflength ) ) 
                end if
             else
                tmpr(i,j)=temp2d(i2,j2)
                tmprmask(i,j)=0
                if (ivar==1.and.k==1) then
                   ssh (i,j)= 0.
                   newdepth (i,j)= 0.
                end if
             end if
             !tmpr(i,j)=temp2d(i2,j2)

             !print *,tmpr(i,j),i2,j2




          end do
          end do

          ! If Mercator unmasks a point for k=1, we should use a nearest neighbour
          ! approach
          if (k==1) then 
              
             hycnearest=0

             do j=1,jdm
             do i=1,idm

                ! if  mercmask is unset and hycom bathy > 1.  , find nearest
                ! neighbour with set mercmask. Set value to that point 
                if (tmprmask(i,j)==0  .and. hycbathy(i,j)>1.) then 

                   ! find nearest point ...
                   ipivh=i
                   jpivh=j
                   call  nearestpoint(hyclon,hyclat,idm,jdm, &
                         hyclon(i,j),hyclat(i,j),ipivh,jpivh,  &
                         a1,a2,a3,a4,hycbathy>1. .and. tmprmask==1,ass,.false.)

                   ! set tmpr to that value
                   if (ass) then

                      ! simply nearest point
                      tmpr(i,j)=tmpr(ipivh,jpivh)
                      if (ivar==1) then
                         ssh (i,j)=ssh (ipivh,jpivh)
                         newdepth (i,j)=newdepth (ipivh,jpivh)
                      end if

                      hycnearest(i,j)=1
                      maxk(i,j)=k
                   else
                      print *,'nearestpoint failed ! '
                      print *,i,j,ipivh,jpivh
                      stop '(nest_offline_mercator)'
                   end if


                elseif (tmprmask(i,j)/=0 .and.hycbathy(i,j)<1. ) then
                   maxk(i,j)=0;
                end if

             end do
             end do

          else if (k>1) then

             ! if  mercmask is unset and hycom bathy > zlev , use old value 
             ! above this point from hycom grid
             do j=1,jdm
             do i=1,idm
                !if (tmprmask(i,j)==0 .and. hycbathy(i,j)>zlev(1)) then
                if (nint(tmprmask(i,j))==0) then
                   tmpr(i,j)=tmprold(i,j)
                else
                   maxk(i,j)=k
                end if
             end do
             end do

          end if


          ! ubavg and vbavg
          if (ivar==3.or.ivar==4) then

             do j=1,jdm
             do i=1,idm
                if (k==maxk(i,j)) then

                   ! z level is midpoint or bottom? Go for bottom for now
                   if (k==1) then
                      dptmp=zlev(1)
                   else
                      dptmp=zlev(k)-zlev(k-1)
                   end if

                   if (ivar==3) then
                      ubavg(i,j)=ubavg(i,j)+tmpr(i,j)*dptmp
                   elseif (ivar==4) then
                      vbavg(i,j)=vbavg(i,j)+tmpr(i,j)*dptmp
                   end if
                end if
             end do
             end do

          end if


          ! The final thing which can happen is that mercator values are defined
          ! at points deeper than hycombathy. For now we neglect this, it is
          ! handled in the vremap routines


          ! Write to temporary file
          call zaiowr(tmpr    ,ip,.false.,amin,amax,99,.false.) !plon
          write(98,104) chfld,k,zlev(k),amin,amax
          !call zaiowr(tmprmask,ip,.false.,amin,amax,99,.false.) !plon

          ! This keeps the old valid values
          where (tmprmask==1.or.hycnearest==1) tmprold=tmpr

       end do
       end do


       do j=1,jdm
       do i=1,idm
       
          if (maxk(i,j)>0)  then
             ubavg(i,j)=ubavg(i,j)/zlev(maxk(i,j))
             vbavg(i,j)=vbavg(i,j)/zlev(maxk(i,j))
          end if
       end do
       end do

       call zaiowr(ubavg   ,ip,.false.,amin,amax,99,.false.) !plon
       write(98,104) 'ubavg   ',0,0,amin,amax

       call zaiowr(vbavg   ,ip,.false.,amin,amax,99,.false.) !plon
       write(98,104) 'vbavg   ',0,0,amin,amax

       call zaiowr(ssh     ,ip,.false.,amin,amax,99,.false.) !plon
       write(98,104) 'ssh     ',0,0,amin,amax

       call zaiowr(newdepth     ,ip,.false.,amin,amax,99,.false.) !plon
       write(98,104) 'newdepth',0,0,amin,amax

       tmpr=maxk
       !print *,tmpr
       call zaiowr(tmpr    ,ip,.false.,amin,amax,99,.false.) !plon
       write(98,104) 'maxk    ',0,0,amin,amax

       ! Diag stuff
       tmpr=hycnearest
       !print *,tmpr
       call zaiowr(tmpr    ,ip,.false.,amin,amax,99,.false.) !plon
       write(98,104) 'nearest ',0,0,amin,amax

       close(98)
       call zaiocl(99)

       ierr=nf90_close(ncid)
       ierr=nf90_close(ncidU)
       ierr=nf90_close(ncidV)

       ! Dump newdepth to a file readable by the nestbat routines
       if (idm > 999 .or. jdm > 999 ) then
          write(tag11,'(i5.5,a,i5.5)')idm,'x',jdm
       else
          write(tag11,'(i3.3,a,i3.3)')idm,'x',jdm
       end if
       allocate(iofld(idm,jdm))
       iofld=newdepth
       ! Close outer boundary
       !iofld(1,:)=0.
       !iofld(:,1)=0.
       !iofld(idm,:)=0.
       !iofld(:,jdm)=0.
       open (unit=10,file='mercatordepths'//trim(tag11)//'.uf',status='replace',form='unformatted')
       write(10) iofld
       close(10)




   104 format(a8," k=",i3," level=",f10.2,"  "," min/max:",2e14.6)



      end program mercator_horizontal_interpolation





         


