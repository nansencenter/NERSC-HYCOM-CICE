module mod_mercatorgrid

integer :: dimlenx, dimleny, dimlenz
real, allocatable, dimension(:,:) :: merclon , merclat ,mercbathy , mercmask  ! T-cell centers
real, allocatable, dimension(:,:) :: merclonu, merclatu, mercmasku ! U-cell
real, allocatable, dimension(:,:) :: merclonv, merclatv, mercmaskv ! V-cell

character(len=*), parameter :: clat  ='gphit', clon ='glamt'
character(len=*), parameter :: clatu ='gphiu', clonu='glamu'
character(len=*), parameter :: clatv ='gphiv', clonv='glamv'
character(len=*), parameter :: cbathy='hdept'
character(len=*), parameter :: cmask ='tmask'
character(len=*), parameter :: cmasku='umask'
character(len=*), parameter :: cmaskv='vmask'
character(len=*), parameter :: bathyfile='ext-mesh_mask.nc'

! Pivot points for hycom
integer, parameter :: numnear=8 !  # points used in gaussian weight
integer, dimension(:,:,:), allocatable :: inear,inearu,inearv,jnear,jnearu,jnearv
real   , dimension(:,:,:), allocatable :: distnear, distnearu, distnearv

contains

       subroutine mercgrid_init(idm,jdm)
       use netcdf
       implicit none

       integer, intent(in) :: idm,jdm

       integer :: ierr
       integer :: k
       integer :: ncid
       integer :: dimidx, dimidy, dimidz
       integer :: tempid,zvarid,lonid,latid, tmpid

       integer, allocatable :: ip(:,:),tmpint(:,:)
       integer, dimension(nf90_max_var_dims) :: dimIDs
       integer :: i2,j2,i,j
       integer :: ndims
       real :: zlev(1), amax, amin 
       logical :: skipip
       real, external :: spherdist

       ! Read Mercator netcdf file
       ierr=NF90_OPEN(bathyfile,NF90_NOWRITE,ncid)

       ! Get dims
       ierr=nf90_inq_dimid(ncid,'x',dimidx)
       ierr=nf90_inquire_dimension(ncid,dimidx,len=dimlenx)
       if (ierr/=NF90_NOERR) stop '(mercgrid:error getting x dim)'
       !print *,dimidx,dimlenx

       ierr=nf90_inq_dimid(ncid,'y',dimidy)
       ierr=nf90_inquire_dimension(ncid,dimidy,len=dimleny)
       if (ierr/=NF90_NOERR) stop '(mercgrid:error getting y dim)'
       !print *,dimidy,dimleny

       ierr=nf90_inq_dimid(ncid,'z',dimidz)
       ierr=nf90_inquire_dimension(ncid,dimidz,len=dimlenz)
       if (ierr/=NF90_NOERR) stop '(mercgrid:error getting z dim)'
       !print *,dimidz,dimlenz
       print *,'Mercator grid size (nx,ny,nz) :',dimlenx, dimleny, dimlenz


       ! Start allocating
       allocate(merclon  (dimlenx,dimleny))
       allocate(merclonu (dimlenx,dimleny))
       allocate(merclonv (dimlenx,dimleny))
       allocate(merclat  (dimlenx,dimleny))
       allocate(merclatu (dimlenx,dimleny))
       allocate(merclatv (dimlenx,dimleny))
       allocate(mercmask (dimlenx,dimleny))
       allocate(mercmasku(dimlenx,dimleny)) ! NB - no z dim
       allocate(mercmaskv(dimlenx,dimleny)) ! NB - no z dim
       allocate(mercbathy(dimlenx,dimleny)) ! NB - no z dim

       allocate(inear (idm,jdm,numnear)) 
       allocate(inearu(idm,jdm,numnear)) 
       allocate(inearv(idm,jdm,numnear)) 

       allocate(jnear (idm,jdm,numnear)) 
       allocate(jnearu(idm,jdm,numnear)) 
       allocate(jnearv(idm,jdm,numnear)) 

       allocate(distnear (idm,jdm,numnear)) 
       allocate(distnearu(idm,jdm,numnear)) 
       allocate(distnearv(idm,jdm,numnear)) 

!!!!!!!!!! Longitudes

       ! Get Mercator lon (T)
       ierr=nf90_inq_varid(ncid,trim(clon),lonid)
       if (ierr==NF90_NOERR) then
          ierr=nf90_get_var(ncid,lonid,merclon)
          if (ierr/=NF90_NOERR) stop '(mercgrid:error getting t lon)'
          print *,'merc lon (T):',minval(merclon),maxval(merclon)
       else
          print *,'Error on reading mercator grid - t lon'
          call exit(1)
       end if

       ! Get Mercator lon (U)
       ierr=nf90_inq_varid(ncid,trim(clonu),lonid)
       if (ierr==NF90_NOERR) then
          ierr=nf90_get_var(ncid,lonid,merclonu)
          if (ierr/=NF90_NOERR) stop '(mercgrid:error getting u lon)'
          print *,'merc lon (U):',minval(merclonu),maxval(merclonu)
       else
          print *,'Error on reading mercator grid - u lon'
          call exit(1)
       end if

       ! Get Mercator lon (V)
       ierr=nf90_inq_varid(ncid,trim(clonv),lonid)
       if (ierr==NF90_NOERR) then
          ierr=nf90_get_var(ncid,lonid,merclonv)
          if (ierr/=NF90_NOERR) stop '(mercgrid:error getting V lon)'
          print *,'merc lon (V):',minval(merclonv),maxval(merclonv)
       else
          print *,'Error on reading mercator grid - V lon'
          call exit(1)
       end if

!!!!!!!!!! Latitudes

       ! Get Mercator lat (T)
       ierr=nf90_inq_varid(ncid,trim(clat),latid)
       if (ierr==NF90_NOERR) then
          ierr=nf90_get_var(ncid,latid,merclat)
          if (ierr/=NF90_NOERR) stop '(mercgrid:error getting t lat)'
          print *,'merc lat (T):',minval(merclat),maxval(merclat)
       else
          print *,'Error on reading mercator grid - t lat'
          call exit(1)
       end if

       ! Get Mercator lat (U)
       ierr=nf90_inq_varid(ncid,trim(clatu),latid)
       if (ierr==NF90_NOERR) then
          ierr=nf90_get_var(ncid,latid,merclatu)
          if (ierr/=NF90_NOERR) stop '(mercgrid:error getting U lat)'
          print *,'merc lat (U):',minval(merclatu),maxval(merclatu)
       else
          print *,'Error on reading mercator grid - U lat'
          call exit(1)
       end if

       ! Get Mercator lat (V)
       ierr=nf90_inq_varid(ncid,trim(clatv),latid)
       if (ierr==NF90_NOERR) then
          ierr=nf90_get_var(ncid,latid,merclatv)
          if (ierr/=NF90_NOERR) stop '(mercgrid:error getting V lat)'
          print *,'merc lat (V):',minval(merclatv),maxval(merclatv)
       else
          print *,'Error on reading mercator grid - V lat'
          call exit(1)
       end if

!!!!!!!!!! The rest - except masks and pivots

       ! Get Mercator bathymetry (T)
       ierr=nf90_inq_varid(ncid,trim(cbathy),tmpid)
       if (ierr==NF90_NOERR) then
          ierr=nf90_get_var(ncid,tmpid,mercbathy)
          if (ierr/=NF90_NOERR) stop '(mercgrid:error getting t bathy)'
          print *,'merc bathy (T):',minval(mercbathy),maxval(mercbathy)
       else
          print *,'Error on reading mercator grid - t bathy'
          call exit(1)
       end if

       ierr=nf90_close(ncid)
       !stop

    end subroutine

    subroutine mercgrid_mask(k)
    use netcdf
    implicit none
    integer, intent(in) :: k ! vertical level

    integer :: ierr, tmpid, ncid
 
    ! Read Mercator netcdf file
    ierr=NF90_OPEN(bathyfile,NF90_NOWRITE,ncid)

    ! Get Mercator mask (T)
    ierr=nf90_inq_varid(ncid,trim(cmask),tmpid)
    if (ierr==NF90_NOERR) then
       ierr=nf90_get_var(ncid,tmpid,mercmask,start=(/1,1,k,1/))
       if (ierr/=NF90_NOERR) stop '(mercgrid_mask:error getting t mask)'
       !print *,'merc mask (T):',minval(mercmask),maxval(mercmask)
    else
       print *,'Error on reading mercator grid - T mask'
       call exit(1)
    end if

    ! Get Mercator mask (U)
    ierr=nf90_inq_varid(ncid,trim(cmasku),tmpid)
    if (ierr==NF90_NOERR) then
       ierr=nf90_get_var(ncid,tmpid,mercmasku,start=(/1,1,k,1/))
       if (ierr/=NF90_NOERR) stop '(mercgrid_mask:error getting u mask)'
       !print *,'merc mask (U):',minval(mercmasku),maxval(mercmasku)
    else
       print *,'Error on reading mercator grid - T mask'
       call exit(1)
    end if

    ! Get Mercator mask (V)
    ierr=nf90_inq_varid(ncid,trim(cmaskv),tmpid)
    if (ierr==NF90_NOERR) then
       ierr=nf90_get_var(ncid,tmpid,mercmaskv,start=(/1,1,k,1/))
       if (ierr/=NF90_NOERR) stop '(mercgrid_mask:error getting v mask)'
       !print *,'merc mask (V):',minval(mercmaskv),maxval(mercmaskv)
    else
       print *,'Error on reading mercator grid - V mask'
       call exit(1)
    end if

    end subroutine mercgrid_mask



    subroutine mercgrid_pivots(merclon,merclat,hyclon,hyclat,inear,jnear,distnear,idm,jdm,usefile,c1)
    use mod_za, only : zaiowr, zaiocl, zaiopf, zaiord
    implicit none
    integer, intent(in) :: idm,jdm
    real, dimension(dimlenx,dimleny), intent(in) :: merclon, merclat ! NB - overrides module vars
    real,    dimension(idm,jdm), intent(in) :: hyclon, hyclat
    real,    dimension(idm,jdm,numnear), intent(out) :: distnear
    integer, dimension(idm,jdm,numnear), intent(out) :: inear,jnear
    character(len=1), intent(in) :: c1
    real, external :: spherdist


    logical, intent(in) :: usefile

    logical ::skipip, match
    integer :: i,j,i2,j2
    real :: dist, mindist,tmpr(idm,jdm),amin,amax
    integer :: ip(idm,jdm)

    real :: mindistnear(numnear)
    integer :: idist, idist2
    integer :: mini(numnear), minj(numnear)

    integer :: k
    logical :: ex

       inquire(exist=ex,file='pivots'//c1//'.a')

       ! Go through hycom grid  - find point closest point on mercator grid
       if (usefile.and.ex) then
          print *,'NB!! Using stored pivots !!! '//c1



          call zaiopf('pivots'//c1//'.a','old',66)

          do k=1,numnear
             call zaiord(tmpr,ip,.false.,amin,amax,66) !plon
             inear(:,:,k)=nint(tmpr)

             call zaiord(tmpr,ip,.false.,amin,amax,66) !plon
             jnear(:,:,k)=nint(tmpr)

             call zaiord(tmpr,ip,.false.,amin,amax,66) !plon
             distnear(:,:,k)=tmpr
          end do
          call zaiocl(66)
       else

          print *,minval(hyclon),maxval(hyclon)
          print *,minval(hyclat),maxval(hyclat)
          print *,minval(merclon),maxval(merclon)
          print *,minval(merclat),maxval(merclat)
          inear=0
          jnear=0
          distnear=0.

          do j=1,jdm
          if (mod(j,jdm/10)==0.or.j==1) print '(a,i4,a,i4,a)',c1//' processed ',j, &
                                        ' of ', jdm,' lines'
          do i=1,idm

             mindist=1e14
             mindistnear=1e14
             do j2=1,dimleny
             do i2=1,dimlenx
                dist=spherdist(hyclon(i,j),hyclat(i,j),merclon(i2,j2),merclat(i2,j2))
                if (dist < mindist) then
                   mindist=dist
                   !print *,i2,j2,mindist
                end if

                if (any(mindistnear>dist)) then
                   match=.false.
                   idist=0
                   do while (.not.match)
                      idist=idist+1
                      if (dist < mindistnear(idist)) match=.true.
                   end do

                   !print *,'idist,dist is ',idist,dist
                   !print *,'idist,dist is ',idist,dist
                   !print *,'mindistnear   is ',mindistnear
                   do idist2=numnear,max(idist,2),-1

                      ! Move stuff down a step
                      if (idist2<numnear) then
                         mini    (idist2)=mini    (idist2-1)
                         minj    (idist2)=minj    (idist2-1)
                         mindistnear(idist2)=mindistnear(idist2-1)
                      end if

                   end do
                   ! This one should be placed here
                   mini    (idist)=i2
                   minj    (idist)=j2
                   mindistnear(idist)=dist
                end if
                      
             end do
             end do
             !print *,'mindist4   is ',mindist4
             !print *,'sum mindist4 is ',sum(mindist4)
             !print *,'mini       is ',mini
             !print *,'minj       is ',minj

             !print *,i,j,ipivtmp,jpivtmp,hyclon(i,j),hyclat(i,j),merclon(ipivtmp,jpivtmp),merclat(ipivtmp,jpivtmp)

             distnear(i,j,:)=mindistnear
             inear   (i,j,:)=mini
             jnear   (i,j,:)=minj
             !print *
             !print *,hyclon(i,j),hyclat(i,j)
             !print *,merclon(mini(1),minj(1)),merclat(mini(1),minj(1))
             !print *,minval(inear(i,j,:)),maxval(inear(i,j,:))
             !print *,minval(jnear(i,j,:)),maxval(jnear(i,j,:))
             !print *,minval(distnear(i,j,:)),maxval(distnear(i,j,:))

          end do
          !print *,minval(inear,inear/=0),maxval(inear)
          !print *,minval(jnear,inear/=0),maxval(jnear)
          !print *,minval(distnear,inear/=0),maxval(distnear)
          end do


          call zaiopf('pivots'//c1//'.a','replace',66)
          do k=1,numnear
             tmpr=real(inear(:,:,k))
             call zaiowr(tmpr,ip,.false.,amin,amax,66,.false.) !plon

             tmpr=real(jnear(:,:,k))
             call zaiowr(tmpr,ip,.false.,amin,amax,66,.false.) !plon

             tmpr=distnear(:,:,k)
             call zaiowr(tmpr,ip,.false.,amin,amax,66,.false.) !plon
          end do




          call zaiocl(66)
       endif

    end subroutine mercgrid_pivots




end module mod_mercatorgrid
