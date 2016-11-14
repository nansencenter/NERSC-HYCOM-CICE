module mod_mercatorgrid

integer :: dimlenx, dimleny, dimlenz
real, allocatable, dimension(:,:) :: merclon , merclat  ! T-cell centers
real, allocatable, dimension(:,:) :: merclonu, merclatv ! U-cell
real, allocatable, dimension(:,:) :: merclonv, merclatv ! V-cell

contains

       subroutine mercgrid(cfld,fld,level)
       use netcdf
       use m_spherdist
       implicit none

       integer, intent(in) :: level
       character(len=*), intent(in) :: cfld
       real, dimension(dimlenx,dimleny) :: fld
       character(len=*), parameter :: bathyfile='ext-mesh_mask.nc'
       character(len=*), parameter :: clon='nav_lon'
       character(len=*), parameter :: clat='nav_lat'
       character(len=*), parameter :: ctmask='tmask'
       character(len=*), parameter :: cbathy='mbathy'


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

       ! Read Mercator netcdf file
       ierr=NF90_OPEN(bathyfile,NF90_NOWRITE,ncid)

       ! Get dims
       ierr=nf90_inq_dimid(ncid,'x',dimidx)
       ierr=nf90_inquire_dimension(ncid,dimidx,len=dimlenx)
       if (ierr/=NF90_NOERR) stop '(mercgrid:error getting x dim)'
       print *,dimidx,dimlenx

       ierr=nf90_inq_dimid(ncid,'y',dimidy)
       ierr=nf90_inquire_dimension(ncid,dimidy,len=dimleny)
       if (ierr/=NF90_NOERR) stop '(mercgrid:error getting y dim)'
       print *,dimidy,dimleny

       ierr=nf90_inq_dimid(ncid,'z',dimidz)
       ierr=nf90_inquire_dimension(ncid,dimidz,len=dimlenz)
       if (ierr/=NF90_NOERR) stop '(mercgrid:error getting z dim)'
       print *,dimidz,dimlenz

       ! Get Mercator lon
       if (trim(cfld)==trim(clon)) then
          ierr=nf90_inq_varid(ncid,trim(clon),lonid)
          ierr=nf90_get_var(ncid,lonid,fld)
          if (ierr/=NF90_NOERR) stop '(mercgrid:error getting v lon)'
          print *,'merc lon:',minval(fld),maxval(fld)
          ierr=nf90_close(ncid)

       ! Get Mercator lat
       elseif (trim(cfld)==trim(clat)) then
          ierr=nf90_inq_varid(ncid,trim(clat),latid)
          ierr=nf90_get_var(ncid,latid,fld)
          if (ierr/=NF90_NOERR) stop '(mercgrid:error getting v lat)'
          print *,'merc lat:',minval(fld),maxval(fld)


       ! Get Mercator bathymetry
       elseif (trim(cfld)==trim(cbathy)) then
          ierr=NF90_OPEN(bathyfile,NF90_NOWRITE,ncid)
          ierr=nf90_inq_varid(ncid,trim(cbathy),tmpid)
          ierr=nf90_get_var(ncid,tmpid,fld)
          ierr=nf90_close(ncid)
          print *,'merc bathy:',minval(fld),maxval(fld)

       ! Get Mercator mask at a level
       elseif (trim(cfld)==trim(ctmask)) then
          allocate(tmpint(dimlenx,dimleny))
          ierr=NF90_OPEN(bathyfile,NF90_NOWRITE,ncid)
          ierr=nf90_inq_varid(ncid,trim(ctmask),tmpid)
          ierr=nf90_get_var(ncid,tmpid,tmpint,start=(/1,1,level,1/)) ! NB varies with z
          ierr=nf90_close(ncid)
          fld=int(tmpint)
          print *,'merc mask:',minval(fld),maxval(fld)

       else
          print *,'Initialized dimension lengths, but thats all'
       end if
    end subroutine


    subroutine mercgrid_pivots(merclon,merclat,hyclon,hyclat,ipiv,jpiv, &
                               i4,j4,dist4,idm,jdm,usefile)
    use m_spherdist
    use mod_za, only : zaiowr, zaiocl, zaiopf, zaiord
    implicit none
    integer, intent(in) :: idm,jdm
    real, dimension(dimlenx,dimleny) :: merclon, merclat
    real, dimension(idm,jdm) :: hyclon, hyclat
    integer, dimension(idm,jdm) :: ipiv,jpiv
    logical, intent(in) :: usefile
    real, dimension(idm,jdm,4) :: dist4
    integer, dimension(idm,jdm,4) :: i4,j4

    logical ::skipip, match
    integer :: i,j,i2,j2
    real :: dist, mindist,tmpr(idm,jdm),amin,amax
    integer :: ipivtmp,jpivtmp,ip(idm,jdm)

    real :: mindist4(4)
    integer :: idist, idist2
    integer :: mini(4), minj(4)

    integer :: k


       ! Go through hycom grid  - find point closest point on mercator grid
       if (usefile) then
          ipiv=1
          jpiv=1
          print *,'NB!! Using stored pivots !!!'
          print *,'NB!! Using stored pivots !!!'
          print *,'NB!! Using stored pivots !!!'
          print *,'NB!! Using stored pivots !!!'
          print *,'NB!! Using stored pivots !!!'
          print *,'NB!! Using stored pivots !!!'
          print *,'NB!! Using stored pivots !!!'
          print *,'NB!! Using stored pivots !!!'
          print *,'NB!! Using stored pivots !!!'
          print *,'NB!! Using stored pivots !!!'
          print *,'NB!! Using stored pivots !!!'
          call zaiopf('pivots.a','old',66)

          call zaiord(tmpr,ip,.false.,amin,amax,66) !plon
          ipiv=int(tmpr)

          call zaiord(tmpr,ip,.false.,amin,amax,66) !plon
          jpiv=int(tmpr)

          do k=1,4
             call zaiord(tmpr,ip,.false.,amin,amax,66) !plon
             i4(:,:,k)=nint(tmpr)

             call zaiord(tmpr,ip,.false.,amin,amax,66) !plon
             j4(:,:,k)=nint(tmpr)

             call zaiord(tmpr,ip,.false.,amin,amax,66) !plon
             dist4(:,:,k)=tmpr
          end do
          call zaiocl(66)
       else
          do j=1,jdm
          print *,j
          do i=1,idm

             mindist=1e14
             mindist4=1e14
             ipiv(i,j)=-1
             jpiv(i,j)=-1
             do j2=1,dimleny
             do i2=1,dimlenx
                dist=spherdist(hyclon(i,j),hyclat(i,j),merclon(i2,j2),merclat(i2,j2))
                if (dist < mindist) then
                   ipivtmp=i2
                   jpivtmp=j2
                   mindist=dist
                   !print *,i2,j2,mindist
                end if

                if (any(mindist4>dist)) then
                   match=.false.
                   idist=0
                   do while (.not.match)
                      idist=idist+1
                      if (dist < mindist4(idist)) match=.true.
                   end do

                   !print *,'idist,dist is ',idist,dist
                   !print *,'mindist4   is ',mindist4
                   do idist2=4,max(idist,2),-1

                      ! Move stuff down a step
                      if (idist2<4) then
                         mini    (idist2)=mini    (idist2-1)
                         minj    (idist2)=minj    (idist2-1)
                         mindist4(idist2)=mindist4(idist2-1)
                      end if

                   end do
                   ! This one should be placed here
                   mini    (idist)=i2
                   minj    (idist)=j2
                   mindist4(idist)=dist
                end if
                      
             end do
             end do
             !print *,'mindist4   is ',mindist4
             !print *,'sum mindist4 is ',sum(mindist4)
             !print *,'mini       is ',mini
             !print *,'minj       is ',minj

             ipiv(i,j)=ipivtmp
             jpiv(i,j)=jpivtmp
             !print *,i,j,ipivtmp,jpivtmp,hyclon(i,j),hyclat(i,j),merclon(ipivtmp,jpivtmp),merclat(ipivtmp,jpivtmp)

             dist4(i,j,:)=mindist4
             i4   (i,j,:)=mini
             j4   (i,j,:)=minj

          end do
          end do
          call zaiopf('pivots.a','replace',66)

          tmpr=real(ipiv)
          call zaiowr(tmpr,ip,.false.,amin,amax,66,.false.) !plon
          tmpr=real(jpiv)
          call zaiowr(tmpr,ip,.false.,amin,amax,66,.false.) !plon

          do k=1,4
             tmpr=real(i4(:,:,k))
             call zaiowr(tmpr,ip,.false.,amin,amax,66,.false.) !plon

             tmpr=real(j4(:,:,k))
             call zaiowr(tmpr,ip,.false.,amin,amax,66,.false.) !plon

             tmpr=dist4(:,:,k)
             call zaiowr(tmpr,ip,.false.,amin,amax,66,.false.) !plon
          end do




          call zaiocl(66)
       endif

    end subroutine mercgrid_pivots




end module mod_mercatorgrid
